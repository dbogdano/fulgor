#pragma once

#include "include/index.hpp"
#include "include/GGCAT.hpp"

namespace fulgor {

struct buffer {
    buffer(uint64_t capacity) : m_capacity(capacity), m_size(0), m_num_sets(0) {
        m_buffer.resize(capacity);
    }

    bool insert(uint32_t const* data, uint32_t size) {
        if (m_size + size + 1 > m_capacity) return false;
        memcpy(m_buffer.data() + m_size, &size, sizeof(uint32_t));
        memcpy(m_buffer.data() + m_size + 1, data, sizeof(uint32_t) * size);
        m_size += size + 1;
        m_num_sets++;
        return true;
    }

    uint32_t size() const { return m_size; }
    uint64_t capacity() const { return m_capacity; }
    uint32_t num_sets() const { return m_num_sets; }
    uint32_t operator[](uint32_t i) const { return m_buffer[i]; }
    uint32_t const* data() const { return m_buffer.data(); }

    void clear() {
        m_size = 0;
        m_num_sets = 0;
    }

private:
    uint64_t m_capacity, m_size, m_num_sets;
    std::vector<uint32_t> m_buffer;
};

template <typename ColorSets>
struct index<ColorSets>::builder {
    builder() {}

    builder(build_configuration const& build_config) : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index already built");

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        {
            essentials::logger("step 1. build colored compacted dBG");
            timer.start();
            m_ccdbg.build(m_build_config);
            m_build_config.num_colors = m_ccdbg.num_colors();
            timer.stop();
            std::cout << "** building the ccdBG took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        std::string input_filename_for_sshash = m_build_config.tmp_dirname + "/" +
                                                util::filename(m_build_config.file_base_name) +
                                                ".sshash.fa";

        {
            essentials::logger("step 2. build m_u2c and m_color_sets");
            timer.start();

            uint64_t num_unitigs = 0;
            uint64_t num_distinct_color_sets = 0;

            typename ColorSets::builder main_builder(m_build_config.num_colors);
            // main_builder.reserve_num_bits(16 * essentials::GB * 8);

            const uint64_t num_threads = m_build_config.num_threads;
            std::vector<typename ColorSets::builder> thread_builders(num_threads,
                                                                     m_build_config.num_colors);

            constexpr uint64_t MAX_BUFFER_SIZE = 1 << 28;
            uint64_t buffer_size = std::min(m_build_config.num_colors * 10000, MAX_BUFFER_SIZE);
            std::vector<std::thread> threads(num_threads);
            std::vector<buffer> thread_buffers(num_threads, buffer_size);

            /* reserve for each build as much space as for the uncompressed buffers */
            for (auto& b : thread_builders) b.reserve_num_bits(buffer_size * 32 / 8);

            assert(thread_buffers[0].capacity() > m_build_config.num_colors);
            uint32_t curr_thread = 0;
            std::atomic<uint32_t> appending_thread = 0;

            auto encode_color_sets_and_append = [&](uint64_t thread_id) {
                buffer const& b = thread_buffers[thread_id];
                thread_builders[thread_id].clear();
                for (uint32_t i = 0, pos = 0; i < b.num_sets(); i++) {
                    uint32_t size = b[pos++];
                    thread_builders[thread_id].encode_color_set(b.data() + pos, size);
                    pos += size;
                }
                while (appending_thread != thread_id) {}
                main_builder.append(thread_builders[thread_id]);
                appending_thread = (appending_thread + 1) % num_threads;
            };

            bits::bit_vector::builder u2c_builder;

            /* write unitigs to fasta file for SSHash */
            std::ofstream out(input_filename_for_sshash.c_str());
            if (!out.is_open()) throw std::runtime_error("cannot open output file");

            m_ccdbg.loop_through_unitigs([&](ggcat::Slice<char> const unitig,
                                             ggcat::Slice<uint32_t> const color_set,
                                             bool same_color_set) {
                assert(curr_thread >= 0);
                assert(curr_thread < num_threads);
                try {
                    if (!same_color_set) {
                        num_distinct_color_sets += 1;
                        if (num_unitigs > 0) u2c_builder.set(num_unitigs - 1, 1);

                        /* fill buffers */
                        if (!thread_buffers[curr_thread].insert(color_set.data, color_set.size)) {
                            threads[curr_thread] =
                                std::thread(encode_color_sets_and_append, curr_thread);
                            const uint32_t next_thread = (curr_thread + 1) % num_threads;
                            if (threads[next_thread].joinable()) threads[next_thread].join();

                            curr_thread = next_thread;

                            thread_buffers[curr_thread].clear();
                            thread_buffers[curr_thread].insert(color_set.data, color_set.size);
                        }
                    }
                    u2c_builder.push_back(0);

                    /*
                        Rewrite unitigs in color-set order.
                        This is *not* the same order in which
                        unitigs are written in the ggcat.fa file.
                    */
                    out << ">\n";
                    out.write(unitig.data, unitig.size);
                    out << '\n';

                    num_unitigs += 1;

                } catch (std::exception const& e) {
                    std::cerr << e.what() << std::endl;
                    exit(1);
                }
            });

            threads[curr_thread] = std::thread(encode_color_sets_and_append, curr_thread);
            for (auto& t : threads) {
                if (t.joinable()) t.join();
            }

            out.close();

            assert(num_unitigs > 0);
            assert(num_unitigs < (uint64_t(1) << 32));

            std::cout << "num_unitigs " << num_unitigs << std::endl;
            std::cout << "num_distinct_color_sets " << num_distinct_color_sets << std::endl;

            main_builder.build(idx.m_color_sets);

            timer.stop();
            std::cout << "** building color sets took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();

            timer.start();
            u2c_builder.set(num_unitigs - 1, 1);
            u2c_builder.build(idx.m_u2c);
            idx.m_u2c_rank1_index.build(idx.m_u2c);
            assert(idx.m_u2c.num_bits() == num_unitigs);
            assert(idx.m_u2c_rank1_index.num_ones() == num_distinct_color_sets);

            std::cout << "m_u2c.num_bits() " << idx.m_u2c.num_bits() << std::endl;
            std::cout << "m_u2c_rank1_index.num_ones() " << idx.m_u2c_rank1_index.num_ones()
                      << std::endl;

            timer.stop();
            std::cout << "** building m_u2c " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 3. build m_k2u");
            timer.start();

            sshash::build_configuration sshash_config;
            sshash_config.k = m_build_config.k;
            sshash_config.m = m_build_config.m;
            sshash_config.canonical = true;
            sshash_config.verbose = m_build_config.verbose;
            sshash_config.tmp_dirname = m_build_config.tmp_dirname;
            sshash_config.num_threads = util::largest_power_of_2(m_build_config.num_threads);
            sshash_config.print();
            idx.m_k2u.build(input_filename_for_sshash, sshash_config);
            try {  // remove unitig file
                std::remove(input_filename_for_sshash.c_str());
            } catch (std::exception const& e) { std::cerr << e.what() << std::endl; }

            timer.stop();
            std::cout << "** building m_k2u took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 4. write filenames");
            timer.start();
            idx.m_filenames.build(m_ccdbg.filenames());
            timer.stop();
            std::cout << "** writing filenames took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        if (m_build_config.check)  //
        {
            essentials::logger("step 5. check correctness...");
            m_ccdbg.loop_through_unitigs(
                [&](ggcat::Slice<char> const unitig,         //
                    ggcat::Slice<uint32_t> const color_set,  //
                    bool /* same_color_set */)               //
                {
                    auto lookup_result = idx.m_k2u.lookup_advanced(unitig.data);
                    const uint64_t unitig_id = lookup_result.contig_id;
                    const uint64_t color_id = idx.u2c(unitig_id);
                    for (uint64_t i = 1; i != unitig.size - idx.m_k2u.k() + 1; ++i) {
                        uint64_t got = idx.m_k2u.lookup_advanced(unitig.data + i).contig_id;
                        if (got != unitig_id) {
                            std::cout << "got unitig_id " << got << " but expected " << unitig_id
                                      << std::endl;
                            return;
                        }
                    }
                    auto fwd_it = idx.m_color_sets.color_set(color_id);
                    const uint64_t size = fwd_it.size();
                    if (size != color_set.size) {
                        std::cout << "got color_set size " << size << " but expected "
                                  << color_set.size << std::endl;
                        return;
                    }
                    for (uint64_t i = 0; i != size; ++i, ++fwd_it) {
                        uint32_t ref = *fwd_it;
                        if (ref != color_set.data[i]) {
                            std::cout << "got ref " << ref << " but expected " << color_set.data[i]
                                      << std::endl;
                            return;
                        }
                    }
                },
                m_build_config.num_threads  //
            );
            essentials::logger("DONE!");
        }
    }

private:
    build_configuration m_build_config;
    GGCAT m_ccdbg;
};

}  // namespace fulgor
