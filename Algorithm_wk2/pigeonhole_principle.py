from kmer_index import *
def approximate_matching_with_index(p, t, n, k):
    ## divide the pattern into n+1 segments
    segment_length = int(round(len(p)/(n+1)))
    ## asssue no repeat
    all_matches = set()
    index = Index(t, k)
    total_hits = 0
    ## test for each segments
    for i in range(n+1):
        start = i * segment_length
        # ensure index within bound
        end = min((i+1) * segment_length, len(p))
        # Call index method to find matches
        matches = queryIndex(p[start:end], t, index)
        hits = len(index.query(p[start:end]))
        total_hits += hits
        # for all found mathes
        for m in matches:
            # check if offset is within the segment
            if m < start or m-start+len(p) > len(t):
                continue
            # m-start = the begining of the pattern
            mismatches = 0
            # check all bases before segment if matching
            for j in range(0,start):
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            # check all bases afet segment if matching
            for j in range(end,len(p)):
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            # if the # of mismatches are less than given n, we add the offset found
            if mismatches <= n:
                all_matches.add(m-start)
    return list(all_matches), total_hits


