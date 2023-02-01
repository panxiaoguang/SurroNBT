using FASTX
using BioSequences
using Fire
## run BWA alignment
## bwa aln ../../DataBase/hg38.fa gRNAs.txt |bwa samse -f gRNA.sam ../../DataBase/hg38.fa - gRNAs.txt
function readFasta(fasta::String)
    sequences = Dict{String,String}()
    reader = open(FASTA.Reader, fasta)
    for record in reader
        sequences[FASTA.identifier(record)] = convert(String, FASTA.sequence(record))
    end
    sequences
end

@main function runBWA(; fasta::String="none", reference::String="none", outfile::String="none")
    gRNAs = readFasta(fasta)
    f = open(outfile, "w")
    open(pipeline(`bwa samse $reference - $fasta`, stdin=pipeline(`bwa aln $reference $fasta`), stderr="log.txt")) do IO
        for line in eachline(IO)
            if !startswith(line, "@")
                name, tag, chrom, start, _... = split(line, "\t")
                gRNA = gRNAs[name]
                if tag == "16"
                    println(f, name, "\t", gRNA, "\t", chrom, "\t", start, "\t", parse(Int64, start) + 22, "\t", "-")
                else
                    println(f, name, "\t", gRNA, "\t", chrom, "\t", start, "\t", parse(Int64, start) + 22, "\t", "+")
                end
            end
        end
    end
    close(f)
end