#include <iostream>
#include <filesystem>

#include "external/sshash/external/gz/zip_stream.cpp"
#include "external/sshash/src/build.cpp"
#include "external/sshash/src/dictionary.cpp"
#include "external/sshash/src/info.cpp"

#include "include/index_types.hpp"
#include "src/index.cpp"
#include "src/color_sets.cpp"
#include "external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"

#include "util.cpp"
#include "build.cpp"
#include "permute.cpp"
#include "pseudoalign.cpp"

int help(char* arg0) {
    std::cout << "== Fulgor: a colored de Bruijn graph index "
                 "========================================"
              << std::endl
              << std::endl;

    std::cout << "Usage: " << arg0 << " <tool> ...\n\n";

    std::cout
        << "Tools:\n"
        << "  build              build an index\n"
        << "  pseudoalign        perform pseudoalignment to an index\n"
        << "  verify             verify that index works correctly with current library version\n"
        << "  stats              print index statistics\n"
        << "  print-filenames    print all reference filenames\n"
        << std::endl;

    std::cout << "Advanced tools:\n"
              << "  permute            permute the reference names of an index\n"
              << "  dump               write unitigs and color sets of an index in text format\n"
              << "  color              build a meta- or a diff- or a meta-diff- index\n"
              << std::endl;

    return 1;
}

int main(int argc, char** argv) {
    if (argc < 2) return help(argv[0]);

    auto tool = std::string(argv[1]);

    /* basic tools */
    if (tool == "build") {
        return build(argc - 1, argv + 1);
    } else if (tool == "pseudoalign") {
        return pseudoalign(argc - 1, argv + 1);
    } else if (tool == "verify") {
        return verify(argc - 1, argv + 1);
    } else if (tool == "stats") {
        return stats(argc - 1, argv + 1);
    } else if (tool == "print-filenames") {
        return print_filenames(argc - 1, argv + 1);
    }

    /* advanced tools */
    else if (tool == "permute") {
        return permute(argc - 1, argv + 1);
    } else if (tool == "dump") {
        return dump(argc - 1, argv + 1);
    } else if (tool == "color") {
        return color(argc - 1, argv + 1);
    }

    std::cout << "Unsupported tool '" << tool << "'.\n" << std::endl;

    return help(argv[0]);
}
