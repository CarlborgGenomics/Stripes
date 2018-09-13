from Bio import bgzf
import sys
import gzip 

def do_something(data, out):
    o = bgzf.open(out, "w")
    with gzip.open(data) as f:
        while True:
            line = f.readline()
            l = line.decode()
            if not l:
                break
            if l[0] == "#":
                print(l.strip(), file=o)
                if "FORMAT" in l:
                    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file=o)
                continue
            lt = l.split()
            # lt[4] can have more than one letter: drop ; can be reference (. or <*> ) or alternative (one letter with or without <*>)
            if lt[4] == "<*>": lt[4] = "."
            elif ",<*>" in lt[4]: lt[4] = lt[4].replace(",<*>", "")
            elif "," in lt[4]: continue # We remove the tri-allelic positions
            # Now that lt[4] is correct we use it for the genotype part
            lt[8] = "GT"
            if lt[4] == ".": lt[9] = "0/0"
            else: lt[9] = "1/1"
            line = "\t".join(lt)
            print(line, file=o)            
    o.close()

do_something(snakemake.input[0], snakemake.output[0])
#do_something(sys.argv[1], sys.argv[2])
