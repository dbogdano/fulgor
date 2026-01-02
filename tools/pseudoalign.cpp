#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include "src/ps_full_intersection.cpp"
#include "src/ps_threshold_union.cpp"

using namespace fulgor;

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
                const bool hybrid_keep_best, const bool use_quality, const uint8_t min_kmer_quality,
                std::ofstream& out_file, std::mutex& iomut, std::mutex& ofile_mut,
                const bool verbose)  //
{
    (void)hybrid_keep_best;  // Used only in threshold-union path
    (void)emit_scores;       // Used only in threshold-union path
    (void)use_quality;       // Used only in threshold-union path
    (void)min_kmer_quality;  // Used only in threshold-union path
    std::vector<uint32_t> colors;  // result of pseudoalignment
    std::vector<uint32_t> scores;  // optional per-color scores
    std::stringstream ss;
    uint64_t buff_size = 0;
    constexpr uint64_t buff_thresh = 50;

    auto rg = rparser.getReadGroup();
    while (rparser.refill(rg)) {
        for (auto const& record : rg) {
            switch (algo) {
                case pseudoalignment_algorithm::FULL_INTERSECTION:
                    index.pseudoalign_full_intersection(record.seq, colors);
                    break;
                case pseudoalignment_algorithm::THRESHOLD_UNION:
                    index.pseudoalign_threshold_union(record.seq, colors, threshold,
                                                      emit_scores ? &scores : nullptr,
                                                      hybrid_keep_best,
                                                      use_quality ? &record.qual : nullptr,
                                                      min_kmer_quality);
                    break;
                default:
                    break;
            }
            buff_size += 1;
            if (!colors.empty()) {
                num_mapped_reads += 1;
                ss << record.name << '\t' << colors.size();
                if (emit_scores && algo == pseudoalignment_algorithm::THRESHOLD_UNION) {
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
            num_reads += 1;
            colors.clear();
            scores.clear();

            if (verbose and num_reads > 0 and num_reads % 1000000 == 0) {
                iomut.lock();
                std::cout << "mapped " << num_reads << " reads" << std::endl;
                iomut.unlock();
            }

            if (buff_size > buff_thresh) {
                std::string outs = ss.str();
                ss.str("");
                ofile_mut.lock();
                out_file.write(outs.data(), outs.size());
                ofile_mut.unlock();
                buff_size = 0;
            }
        }
    }

    // dump anything left in the buffer
    if (buff_size > 0) {
        std::string outs = ss.str();
        ss.str("");
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
                const bool hybrid_keep_best, const bool use_quality, const uint8_t min_kmer_quality,
                const bool verbose) {
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

    for (uint64_t i = 1; i != num_threads; ++i) {
        workers.push_back(std::thread([&index, &rparser, &num_reads, &num_mapped_reads, ps_alg,
                                       threshold, emit_scores, hybrid_keep_best, use_quality,
                                       min_kmer_quality, &out_file, &iomut, &ofile_mut,
                                       verbose]() {
            pseudoalign(index, rparser, num_reads, num_mapped_reads, ps_alg, threshold,
                        emit_scores, hybrid_keep_best, use_quality, min_kmer_quality, out_file,
                        iomut, ofile_mut, verbose);
        }));
    }

    for (auto& w : workers) w.join();
    rparser.stop();

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
    parser.add("hybrid_all_hits",
               "For HYBRID indexes, return all colors with scores >= threshold (default: best only)",
               "--hybrid-all-hits", false, true);
    parser.add("use_quality", "Weight k-mers by Phred quality scores (threshold-union only).",
               "--use-quality", false, true);
    parser.add("min_kmer_quality",
               "Minimum average k-mer quality (0-40, default: 0=no filtering). Only used with "
               "--use-quality.",
               "--min-kmer-quality", false);
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
    bool hybrid_keep_best = !parser.get<bool>("hybrid_all_hits");
    bool use_quality = parser.get<bool>("use_quality");
    uint8_t min_kmer_quality = 0;
    if (parser.parsed("min_kmer_quality")) {
        uint32_t val = parser.get<uint32_t>("min_kmer_quality");
        if (val > 40) {
            std::cerr << "min_kmer_quality must be in range [0, 40]" << std::endl;
            return 1;
        }
        min_kmer_quality = static_cast<uint8_t>(val);
    }
    if (min_kmer_quality > 0 && !use_quality) {
        std::cerr << "min_kmer_quality requires --use-quality flag" << std::endl;
        return 1;
    }

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
                                                         ps_alg, emit_scores, hybrid_keep_best,
                                                         use_quality, min_kmer_quality, verbose);
    } else if (sshash::util::ends_with(index_filename,
                                       constants::meta_colored_fulgor_filename_extension)) {
        return pseudoalign<meta_index_type>(index_filename, query_filename, output_filename,
                                            num_threads, threshold, ps_alg, emit_scores,
                                            hybrid_keep_best, use_quality, min_kmer_quality,
                                            verbose);
    } else if (sshash::util::ends_with(index_filename,
                                       constants::diff_colored_fulgor_filename_extension)) {
        return pseudoalign<differential_index_type>(index_filename, query_filename, output_filename,
                                                    num_threads, threshold, ps_alg, emit_scores,
                                                    hybrid_keep_best, use_quality,
                                                    min_kmer_quality, verbose);
    } else if (sshash::util::ends_with(index_filename, constants::fulgor_filename_extension)) {
        return pseudoalign<index_type>(index_filename, query_filename, output_filename, num_threads,
                                       threshold, ps_alg, emit_scores, hybrid_keep_best,
                                       use_quality, min_kmer_quality, verbose);
    }

    std::cerr << "Wrong index filename supplied." << std::endl;

    return 1;
}
