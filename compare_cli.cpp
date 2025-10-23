#include "comparison.h"
#include "rotation_tree.h"
#include "tree_generators.h"

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <string_view>

struct RunnerConfig {
    std::string case_type = "comb";
    int n = 5;
    int count = 1;
    long long seed = 12345;
    double time_limit = 5.0;
    std::size_t visited_cap = 5'000'000;
    std::size_t queue_cap = 5'000'000;
    std::string program = "test_asan";
    bool fallback_bidir = false;
    std::size_t bidir_cap = 5'000'000;
    bool prefer_bidir = false;
};

static void printUsage(const char* argv0) {
    std::cerr << "Usage: " << argv0 << " [--case comb|random] [--n N] [--count C]"
              << " [--seed S] [--time-limit T] [--visited-cap V] [--queue-cap Q]"
              << " [--program NAME] [--fallback-bidir] [--bidir-cap CAP] [--prefer-bidir]\n";
}

static bool parseArguments(int argc, char** argv, RunnerConfig& cfg) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto consume = [&](std::string_view name) -> std::string {
            if (arg == name) {
                if (i + 1 >= argc) {
                    throw std::invalid_argument(std::string(name) + " requires a value");
                }
                return std::string(argv[++i]);
            }
            auto pos = arg.find('=');
            if (pos != std::string::npos && arg.substr(0, pos) == name) {
                return arg.substr(pos + 1);
            }
            return {};
        };

        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            std::exit(0);
        }

        if (auto value = consume("--case"); !value.empty()) {
            cfg.case_type = value;
            continue;
        }
        if (auto value = consume("--n"); !value.empty()) {
            cfg.n = std::stoi(value);
            continue;
        }
        if (auto value = consume("--count"); !value.empty()) {
            cfg.count = std::stoi(value);
            continue;
        }
        if (auto value = consume("--seed"); !value.empty()) {
            cfg.seed = std::stoll(value);
            continue;
        }
        if (auto value = consume("--time-limit"); !value.empty()) {
            cfg.time_limit = std::stod(value);
            continue;
        }
        if (auto value = consume("--visited-cap"); !value.empty()) {
            cfg.visited_cap = static_cast<std::size_t>(std::stoull(value));
            continue;
        }
        if (auto value = consume("--queue-cap"); !value.empty()) {
            cfg.queue_cap = static_cast<std::size_t>(std::stoull(value));
            continue;
        }
        if (auto value = consume("--program"); !value.empty()) {
            cfg.program = value;
            continue;
        }
        if (arg == "--fallback-bidir") {
            cfg.fallback_bidir = true;
            continue;
        }
        if (arg == "--prefer-bidir") {
            cfg.prefer_bidir = true;
            continue;
        }
        if (auto value = consume("--bidir-cap"); !value.empty()) {
            cfg.bidir_cap = static_cast<std::size_t>(std::stoull(value));
            continue;
        }

        std::cerr << "Unknown argument: " << arg << "\n";
        printUsage(argv[0]);
        return false;
    }
    return true;
}

int main(int argc, char** argv) {
    RunnerConfig cfg;
    try {
        if (!parseArguments(argc, argv, cfg)) return 1;
    } catch (const std::exception& ex) {
        std::cerr << "Argument error: " << ex.what() << "\n";
        printUsage(argv[0]);
        return 1;
    }

    if (cfg.n <= 0) {
        std::cerr << "n must be positive\n";
        return 1;
    }
    if (cfg.count <= 0) {
        std::cerr << "count must be positive\n";
        return 1;
    }

    for (int idx = 0; idx < cfg.count; ++idx) {
        VectorRangeTreeMap start, target;
        long long pair_seed = cfg.seed + idx;

        if (cfg.case_type == "comb") {
            Traversals leftComb  = makeCombTraversals(cfg.n, /*rightComb=*/false);
            Traversals rightComb = makeCombTraversals(cfg.n, /*rightComb=*/true);
            start.build(leftComb.preorder, leftComb.inorder);
            target.build(rightComb.preorder, rightComb.inorder);
        } else if (cfg.case_type == "random") {
            std::mt19937 rngA(static_cast<std::uint32_t>(pair_seed * 2 + 0));
            std::mt19937 rngB(static_cast<std::uint32_t>(pair_seed * 2 + 1));
            Traversals tA = makeRandomTraversals(cfg.n, rngA);
            Traversals tB = makeRandomTraversals(cfg.n, rngB);
            start.build(tA.preorder, tA.inorder);
            target.build(tB.preorder, tB.inorder);
        } else {
            std::cerr << "Unsupported case type: " << cfg.case_type << "\n";
            return 1;
        }

        ComparisonOptions opts;
        opts.program = cfg.program;
        opts.case_type = cfg.case_type;
        opts.seed = (cfg.case_type == "random") ? pair_seed : -1;
        opts.time_limit_sec = cfg.time_limit;
        opts.visited_cap = cfg.visited_cap;
        opts.queue_cap = cfg.queue_cap;
        opts.use_bidir_on_timeout = cfg.fallback_bidir;
        opts.bidir_state_cap = cfg.bidir_cap;
        opts.prefer_bidir = cfg.prefer_bidir;

        ComparisonPair pair = runBidirectional(start, target, opts);
        if (!pair.distance_agrees &&
            pair.forward.distance >= 0 && pair.reverse.distance >= 0) {
            std::cerr << "[WARN] distance mismatch for seed " << opts.seed
                      << " n=" << cfg.n << "\n";
        }

        std::cout << comparisonRowToJson(pair.forward) << "\n";
        std::cout << comparisonRowToJson(pair.reverse) << "\n";
    }

    return 0;

    
}
