import sys
import gzip

def process_line(line, populations, anc_index, ibs_ind1_index,prefix):
    fields = line.split('\t')
    anc = fields[anc_index]
    if anc == "0,0:0":
        return None  # Skip line as per requirement
    derived_index = 0 if anc.startswith("0,1") else 1  # if the anc is 1,0 then the derived is index 1 otherwise, the anc is 0,1 and the derived is index 0
    results = {}
    # Process IBS_ind1 separately for base counts
    targetfields = fields[ibs_ind1_index].split(':')[0].split(',')
    
    allele_counts = int(targetfields[derived_index])
    total_counts = int(targetfields[0]) + int(targetfields[1])
    results[prefix] = (allele_counts, total_counts)
    # Process other populations for frequencies
    for pop in populations:
        if pop == prefix:  # Skip IBS_ind1 for frequency calculation
            continue
        targetfields = fields[populations[pop]].split(':')[0].split(',')
        allele_counts = int(targetfields[derived_index])
        total_counts = int(targetfields[0]) + int(targetfields[1])
        if total_counts > 0:  # Avoid division by zero
            freq = float(allele_counts) / total_counts
        else:
            freq = 0.0
        results[pop] = freq
    return results

def main(prefix):
    header = sys.stdin.readline().strip().split('\t')
    populations = {pop: index for index, pop in enumerate(header) if pop not in ['#chr', 'coord', 'REF,ALT', 'root', 'anc']}
    anc_index = header.index('anc')
    ibs_ind1_index = populations[prefix]

    # Open file pointers
    out_files = {pop: gzip.open(f'{prefix}_freq_{pop}.gz', 'wt') for pop in populations if pop != prefix}
    out_files[prefix] = gzip.open(f'{prefix}_basecount_{prefix}.gz', 'wt')

    for line in sys.stdin:
        if line.strip() and not line.startswith('#'):
            results = process_line(line, populations, anc_index, ibs_ind1_index,prefix)
            
            if results:
                #print(line)
                #print(results)
                # Write base counts for IBS_ind1
                allele_counts, total_counts = results[prefix]
                #print((f'{allele_counts} {total_counts}\n'));
                out_files[prefix].write(f'{allele_counts} {total_counts}\n')
                # Write frequencies for other populations
                for pop in results:
                    if pop != prefix:  # Skip IBS_ind1 for frequency file
                        out_files[pop].write(f'{results[pop]}\n')

    # Close file pointers
    for file in out_files.values():
        file.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python script.py <prefix>")
        sys.exit(1)
    main(sys.argv[1])
