from kmer_index import *
def approximate_matching_with_index(p, t, n, k):
    segment_length = int(round(len(p)/(n+1)))
    all_matches = set()
    index = Index(t, k)
    total_hits = 0
    for i in range(n+1):
        start = i * segment_length
        end = min((i+1) * segment_length, len(p))
        matches = queryIndex(p[start:end], t, index)
        hits = len(index.query(p[start:end]))
        total_hits += hits

        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue

            mismatches = 0
            for j in range(0,start):
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            for j in range(end,len(p)):
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                all_matches.add(m-start)
    return list(all_matches), total_hits


