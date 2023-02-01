using Fire
using FASTX
## read in the genomic coordinates
function readCoordinates(filename::String)
    coordinates = Dict{String,Vector{NamedTuple{(:start, :stop, :strand),Tuple{Int64,Int64,String}}}}()
    for line in eachline(filename)
        chrom, start, stop, strand = split(line, '\t')
        if haskey(coordinates, chrom) == false
            coordinates[chrom] = Vector{NamedTuple{(:start, :stop, :strand),Tuple{Int64,Int64,String}}}()
        end
        push!(coordinates[chrom], (start=parse(Int64, start), stop=parse(Int64, stop), strand=strand))
    end
    coordinates
end

## get expanded sequences from the genomic coordinates
function getSequence(coord::Dict{String,Vector{NamedTuple{(:start, :stop, :strand),Tuple{Int64,Int64,String}}}}, fastafile::String, leftFlank::Int64, rightFlank::Int64)
    readers = open(FASTA.Reader, fastafile, index=string(fastafile, ".fai"))
    seqs = Dict{String,String}()
    start = 0
    stop = 0
    expand = dna"NNN"
    for (chrom, regions) in coord
        totalChromSequence = FASTA.sequence(readers[chrom])
        for region in regions
            if region.strand == "+"
                start = region.start - leftFlank
                stop = region.stop + rightFlank
                expand = totalChromSequence[start:stop]
            else
                start = region.start - rightFlank
                stop = region.stop + leftFlank
                expand = reverse_complement(totalChromSequence[start:stop])
            end
            if start < 1 || stop > FASTA.seqlen(readers[chrom])
                println("skip $chrom:$(region.start)-$(region.stop) beacuse it is out of range!")
            else
                seqs["$chrom:$(region.start)-$(region.stop)"] = convert(String, expand)
            end
        end
    end
    seqs
end

@main function makeOntargetLibrary(; input::String="none", output::String="none", reference::String="none", leftFlank=10, rightFlank=4, recog="ACCACGTCTCACACC", scaffold="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT", bind="GTTTGGAGACGACGG")
    f = open(output, "w")
    coord = readCoordinates(input)
    seqs = getSequence(coord, reference, leftFlank, rightFlank)
    for (name, seq) in seqs
        gRNA = seq[11:30]
        chipSequence = string(recog, "g", gRNA, scaffold, seq, bind)
        println(f, name, "\t", chipSequence)
    end
    close(f)
end


