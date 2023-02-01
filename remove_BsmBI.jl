using BioSequences
using Fire
# remove BsnBl from a sequence

@main function removeBsmBl(; fileName::String="none", outFileName::String="none")
    f = open(outFileName, "w")
    BsmBI = ExactSearchQuery(dna"GAGACG")
    BsmBI_REV = ExactSearchQuery(dna"CGTCTC")
    for line in eachline(fileName)
        (id, seq) = split(line, "\t")
        seq = LongDNASeq(seq)
        surroSeq = seq[1:6]
        if occursin(BsmBI, seq) || occursin(BsmBI_REV, seq)
            continue
        else
            println(f, id, "\t", seq)
        end
    end
    close(f)
end
