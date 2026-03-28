#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <map>
#include <mutex>

#include "src/ps_full_intersection.cpp"
#include "src/ps_threshold_union.cpp"

using namespace fulgor;

struct OrderedOutputBuffer {
    std::map<uint64_t, std::string> buffer;
    uint64_t next_write_idx = 0;
    std::mutex mtx;
    std::ofstream* out_file = nullptr;

    void insert(uint64_t idx, const std::string& data) {
        std::lock_guard<std::mutex> lock(mtx);
        buffer[idx] = data;
        flush_ready();
    }

    void flush_ready() {
        // Must be called with lock held
        while (buffer.count(next_write_idx)) {
            out_file->write(buffer[next_write_idx].data(), buffer[next_write_idx].size());
            buffer.erase(next_write_idx);
            next_write_idx++;
        }
    }

    void flush_all() {
        std::lock_guard<std::mutex> lock(mtx);
        for (uint64_t i = next_write_idx; i <= next_write_idx + buffer.size(); ++i) {
            if (buffer.count(i)) {
                out_file->write(buffer[i].data(), buffer[i].size());
            }
        }
        buffer.clear();
    }
};

enum class pseudoalignment_algorithm : uint8_t { FULL_INTERSECTION, THRESHOLD_UNION };

std::string to_string(pseudoalignment_algorithm algo, double threshold) {
    std::string o;
    switch (algo) {
        case pseudoalignment_algorithm::FULL_INTERSECTION:
            o = "full-intersection";
            break;
        case pseudoalignment_algorithm::THRESHOLD_UNION:
            o = "threshold-union (threshold = " + std::to_string(threshold) + ")";
            break;
    }
    return o;
}

template <typename FulgorIndex>
int pseudoalign(FulgorIndex const& index, fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser,
                std::atomic<uint64_t>& num_reads, std::atomic<uint64_t>& num_mapped_reads,
                pseudoalignment_algorithm algo, const double threshold, const bool emit_scores,
                std::ofstream& out_file, std::mutex& ofile_mut,
                OrderedOutputBuffer* ordered_output_buffer, const bool ordered_output,
                std::mutex& iomut,
                const bool verbose)  //
{
    std::vector<uint32_t> colors;  // result of pseudoalignment
    std::vector<uint32_t> scores;  // optional per-color scores
    std::stringstream out_ss;
    uint64_t buff_size = 0;
    constexpr uint64_t buff_thresh = 50;
    uint64_t read_idx = 0;

    auto rg = rparser.getReadGroup();
    while (rparser.refill(rg)) {
        for (auto const& record : rg) {
            read_idx = num_reads++;
            
            switch (algo) {
                case pseudoalignment_algorithm::FULL_INTERSECTION:
                    index.pseudoalign_full_intersection(record.seq, colors);
                    break;
                case pseudoalignment_algorithm::THRESHOLD_UNION:
                    index.pseudoalign_threshold_union(record.seq, colors, threshold,
                                                      emit_scores ? &scores : nullptr);
                    break;
                default:
                    break;
            }
            
            std::stringstream ss;
            if (!colors.empty()) {
                num_mapped_reads += 1;
                ss << record.name << '\t' << colors.size();
                if (emit_scores) {
                    assert(colors.size() == scores.size());
                    for (uint64_t idx = 0; idx != colors.size(); ++idx) {
                        ss << "\t" << colors[idx] << ":" << scores[idx];
                    }
                } else {
                    for (auto c : colors) { ss << "\t" << c; }
                }
                ss << '\n';
            } else {
                ss << record.name << "\t0\n";
            }

            if (ordered_output) {
                ordered_output_buffer->insert(read_idx, ss.str());
            } else {
                out_ss << ss.str();
                buff_size += 1;
                if (buff_size > buff_thresh) {
                    std::string outs = out_ss.str();
                    out_ss.str("");
                    ofile_mut.lock();
                    out_file.write(outs.data(), outs.size());
                    ofile_mut.unlock();
                    buff_size = 0;
                }
            }
            
            colors.clear();
            scores.clear();

            if (verbose && read_idx > 0 && read_idx % 1000000 == 0) {
                iomut.lock();
                std::cout << "mapped " << read_idx << " reads" << std::endl;
                iomut.unlock();
            }
        }
    }

    if (!ordered_output and buff_size > 0) {
        std::string outs = out_ss.str();
        out_ss.str("");
        ofile_mut.lock();
        out_file.write(outs.data(), outs.size());
        ofile_mut.unlock();
        buff_size = 0;
    }

    return 0;
}

template <typename FulgorIndex>
int pseudoalign(std::string const& index_filename, std::string const& query_filename,
                std::string const& output_filename, uint64_t num_threads, double threshold,
                pseudoalignment_algorithm ps_alg, const bool emit_scores,
                const bool ordered_output, const bool verbose) {
    FulgorIndex index;
    if (verbose) essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    if (verbose) essentials::logger("DONE");

    std::cerr << "query mode : " << to_string(ps_alg, threshold) << "\n";

    std::ifstream is(query_filename.c_str());
    if (!is.good()) {
        std::cerr << "error in opening the file '" + query_filename + "'" << std::endl;
        return 1;
    }

    if (verbose) essentials::logger("performing queries from file '" + query_filename + "'...");
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::milliseconds> t;
    t.start();

    std::atomic<uint64_t> num_mapped_reads{0};
    std::atomic<uint64_t> num_reads{0};

    auto query_filenames = std::vector<std::string>({query_filename});
    assert(num_threads >= 2);
    fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(query_filenames, num_threads,
                                                             num_threads - 1);

    rparser.start();
    std::vector<std::thread> workers;
    workers.reserve(num_threads);
    std::mutex iomut;
    std::mutex ofile_mut;

    std::ofstream out_file;
    out_file.open(output_filename, std::ios::out | std::ios::trunc);
    if (!out_file) {
        std::cerr << "could not open output file " + output_filename << std::endl;
        return 1;
    }

    OrderedOutputBuffer output_buffer;
    if (ordered_output) output_buffer.out_file = &out_file;

    for (uint64_t i = 1; i != num_threads; ++i) {
        workers.push_back(std::thread([&index, &rparser, &num_reads, &num_mapped_reads, ps_alg,
                                       threshold, emit_scores, &out_file, &ofile_mut,
                                       &output_buffer, ordered_output, &iomut, verbose]() {
            pseudoalign(index, rparser, num_reads, num_mapped_reads, ps_alg, threshold,
                        emit_scores, out_file, ofile_mut, &output_buffer, ordered_output, iomut,
                        verbose);
        }));
    }

    for (auto& w : workers) w.join();
    rparser.stop();

    // Ensure all buffered output is written
    if (ordered_output) output_buffer.flush_all();

    out_file.close();

    t.stop();
    if (verbose) essentials::logger("DONE");

    if (verbose) {
        std::cout << "mapped " << num_reads << " reads" << std::endl;
        std::cout << "elapsed = " << t.elapsed() << " millisec / ";
        std::cout << t.elapsed() / 1000 << " sec / ";
        std::cout << t.elapsed() / 1000 / 60 << " min / ";
        std::cout << (t.elapsed() * 1000) / num_reads << " musec/read" << std::endl;
        std::cout << "num_mapped_reads " << num_mapped_reads << "/" << num_reads << " ("
                  << (num_mapped_reads * 100.0) / num_reads << "%)" << std::endl;
    }

    return 0;
}

int pseudoalign(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);

    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    parser.add("query_filename", "Query filename in FASTA/FASTQ format (optionally gzipped).", "-q",
               true);
    parser.add("output_filename",
               "File where output will be written. You can specify \"/dev/stdout\" to write "
               "output to stdout. In this case, it is also recommended to use the --verbose flag "
               "to avoid printing status messages to stdout.",
               "-o", true);
    parser.add("num_threads", "Number of threads (default is 1).", "-t", false);
    parser.add("verbose", "Verbose output during query (default is false).", "--verbose", false, true);
    parser.add("threshold",
               "Threshold for threshold_union algorithm. It must be a float in (0.0,1.0].", "-r",
               false);
    parser.add("emit_scores", "Emit color:score pairs for threshold-union output.",
               "--emit-scores", false, true);
    parser.add("ordered_output",
               "Preserve input order in output when using multiple threads (default is false).",
               "--ordered-output", false, true);
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");
    auto query_filename = parser.get<std::string>("query_filename");
    auto output_filename = parser.get<std::string>("output_filename");

    uint64_t num_threads = 1;
    if (parser.parsed("num_threads")) num_threads = parser.get<uint64_t>("num_threads");
    if (num_threads == 1) {
        num_threads += 1;
        std::cerr
            << "1 thread was specified, but an additional thread will be allocated for parsing"
            << std::endl;
    }

    double threshold = constants::invalid_threshold;
    if (parser.parsed("threshold")) threshold = parser.get<double>("threshold");
    if (threshold == 0.0 or threshold > 1.0) {
        std::cerr << "threshold must be a float in (0.0,1.0]" << std::endl;
        return 1;
    }

    bool emit_scores = parser.get<bool>("emit_scores");
    bool ordered_output = parser.get<bool>("ordered_output");

    auto ps_alg = pseudoalignment_algorithm::FULL_INTERSECTION;
    if (threshold != constants::invalid_threshold) {
        ps_alg = pseudoalignment_algorithm::THRESHOLD_UNION;
    }

    bool verbose = parser.get<bool>("verbose");
    if (verbose) util::print_cmd(argc, argv);

    if (sshash::util::ends_with(index_filename,
                                constants::meta_diff_colored_fulgor_filename_extension)) {
        return pseudoalign<meta_differential_index_type>(index_filename, query_filename,
                                                         output_filename, num_threads, threshold,
                                                         ps_alg, emit_scores, ordered_output,
                                                         verbose);
    } else if (sshash::util::ends_with(index_filename,
                                       constants::meta_colored_fulgor_filename_extension)) {
        return pseudoalign<meta_index_type>(index_filename, query_filename, output_filename,
                                            num_threads, threshold, ps_alg, emit_scores,
                                            ordered_output, verbose);
    } else if (sshash::util::ends_with(index_filename,
                                       constants::diff_colored_fulgor_filename_extension)) {
        return pseudoalign<differential_index_type>(index_filename, query_filename, output_filename,
                                                    num_threads, threshold, ps_alg, emit_scores,
                                                    ordered_output, verbose);
    } else if (sshash::util::ends_with(index_filename, constants::fulgor_filename_extension)) {
        return pseudoalign<index_type>(index_filename, query_filename, output_filename, num_threads,
                                       threshold, ps_alg, emit_scores, ordered_output, verbose);
    }

    std::cerr << "Wrong index filename supplied." << std::endl;

    return 1;
}
