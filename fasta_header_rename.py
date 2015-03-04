import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, 'r') as fh_in:
    with open(outfile, 'w') as fh_out:
        result = {}
        count = 1
        for line in fh_in:
            line = line.strip()
            if line.startswith(">"):
                line = line[1:]
                result[line] = "lincRNA_" + str(count)
                count = count + 1
                header = ">" + str(result[line])
                fh_out.write(header)
                fh_out.write("\n")
            else:
                fh_out.write(line)
                fh_out.write("\n")

