// Builds JSON records describing a single solver run. The comparison harness
// (`run_compare.py`) uses these helpers to log distances, timings, queue
// stats, and solver metadata when cross-checking the C++ and Java binaries.
// The module also owns the bidirectional wrapper so we can replay timeouts
// with the hashed meet-in-the-middle search before emitting the final row.

#include "comparison.h"
#include <chrono>
#include <cstdio>
#include <climits>
#include <iomanip>
#include <sstream>
#include <utility>

namespace {

std::string escapeJson(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 16);
    for (unsigned char c : s) {
        switch (c) {
        case '\"': out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\b': out += "\\b"; break;
        case '\f': out += "\\f"; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default:
            if (c < 0x20) {
                char buf[7];
                std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned int>(c));
                out += buf;
            } else {
                out.push_back(static_cast<char>(c));
            }
        }
    }
    return out;
}

} // namespace

ComparisonRow runComparison(const VectorRangeTreeMap& start,
                            const VectorRangeTreeMap& target,
                            const ComparisonOptions& opts)
{
    ComparisonRow row{};
    row.program   = opts.program;
    row.case_type = opts.case_type;
    row.direction = opts.direction;
    row.seed      = opts.seed;
    row.n         = static_cast<int>(start.original_nodes.size());
    row.tree_a    = canonicalTraversalString(start);
    row.tree_b    = canonicalTraversalString(target);
    row.solver    = "bfs";

    if (TreesEqual(start, target)) {
        row.distance = 0;
        row.time_ms = 0.0;
        row.expanded = 0;
        row.enqueued = 0;
        row.visited = start.original_nodes.size();
        row.max_queue = 0;
        row.duplicates = 0;
        row.status = "ok";
        return row;
    }

    if (opts.prefer_bidir) {
        auto fb_start = std::chrono::steady_clock::now();
        int dist = BiBFSSearchHashed(start, target, opts.bidir_state_cap);
        auto fb_end = std::chrono::steady_clock::now();
        row.time_ms = std::chrono::duration<double, std::milli>(fb_end - fb_start).count();
        row.solver = "bidir";
        if (dist != INT_MAX) {
            row.distance = dist;
            row.status = "ok";
        } else {
            row.distance = -1;
            row.status = "fallback_fail";
        }
        return row;
    }

    BFSRun run = BFSSearchCapped(start, target,
                                 opts.time_limit_sec,
                                 opts.visited_cap,
                                 opts.queue_cap);

    row.distance   = (run.dist == INT_MAX) ? -1 : run.dist;
    row.time_ms    = run.seconds * 1000.0;
    row.expanded   = run.expanded;
    row.enqueued   = run.enqueued;
    row.visited    = run.visited;
    row.max_queue  = run.max_queue;
    row.duplicates = run.duplicates;

    if (run.timeout) {
        row.status = "timeout";
    } else if (run.cap_hit) {
        row.status = "cap";
    } else if (run.dist == INT_MAX) {
        row.status = "error:dist_unavailable";
    } else {
        row.status = "ok";
    }

    const bool need_fallback = opts.use_bidir_on_timeout &&
        (row.status == "timeout" || row.status == "cap" || row.distance < 0);
    if (need_fallback) {
        auto fb_start = std::chrono::steady_clock::now();
        int fb_dist = BiBFSSearchHashed(start, target, opts.bidir_state_cap);
        auto fb_end = std::chrono::steady_clock::now();
        row.time_ms = std::chrono::duration<double, std::milli>(fb_end - fb_start).count();
        row.solver = "bidir";
        if (fb_dist != INT_MAX) {
            row.distance = fb_dist;
            row.status = "fallback_ok";
        } else if (row.status != "timeout" && row.status != "cap") {
            row.status = "fallback_fail";
        }
    }

    return row;
}

std::string comparisonRowToJson(const ComparisonRow& row) {
    std::ostringstream out;
    out.setf(std::ios::fixed, std::ios::floatfield);
    out << std::setprecision(3);

    auto toULL = [](std::size_t v) {
        return static_cast<unsigned long long>(v);
    };

    out << "{\"program\":\"" << escapeJson(row.program) << '"'
        << ",\"case_type\":\"" << escapeJson(row.case_type) << '"'
        << ",\"n\":" << row.n
        << ",\"seed\":" << row.seed
        << ",\"direction\":\"" << escapeJson(row.direction) << '"'
        << ",\"distance\":" << row.distance
        << ",\"time_ms\":" << row.time_ms
        << ",\"expanded\":" << toULL(row.expanded)
        << ",\"enqueued\":" << toULL(row.enqueued)
        << ",\"visited\":" << toULL(row.visited)
        << ",\"max_queue\":" << toULL(row.max_queue)
        << ",\"duplicates\":" << toULL(row.duplicates)
        << ",\"status\":\"" << escapeJson(row.status) << '"'
        << ",\"solver\":\"" << escapeJson(row.solver) << '"'
        << ",\"tree_a\":\"" << escapeJson(row.tree_a) << '"'
        << ",\"tree_b\":\"" << escapeJson(row.tree_b) << "\"}";
    return out.str();
}

ComparisonPair runBidirectional(const VectorRangeTreeMap& a,
                                const VectorRangeTreeMap& b,
                                const ComparisonOptions& opts) {
    ComparisonOptions forward_opts = opts;
    forward_opts.direction = "a->b";
    ComparisonRow forward = runComparison(a, b, forward_opts);

    ComparisonOptions reverse_opts = opts;
    reverse_opts.direction = "b->a";
    ComparisonRow reverse = runComparison(b, a, reverse_opts);

    auto solved = [](const ComparisonRow& row) {
        return row.distance >= 0 &&
               (row.status == "ok" || row.status == "fallback_ok");
    };

    bool agree = false;
    if (solved(forward) && solved(reverse)) {
        agree = (forward.distance == reverse.distance);
    }

    return ComparisonPair{std::move(forward), std::move(reverse), agree};
}
