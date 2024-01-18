using BioSequences
using Fire
using FASTX
function makeOneMismatch(gRNA::LongDNASeq)
    OneMismatches = Vector{LongDNASeq}()
    for (idx, base) in enumerate(gRNA)
        for b in [DNA_A, DNA_C, DNA_G, DNA_T]
            if b != base
                newgRNA = copy(gRNA)
                newgRNA[idx] = b
                push!(OneMismatches, newgRNA)
            end
        end
    end
    OneMismatches
end

function makeTwoMismatch(gRNA::LongDNASeq)
    TwoMismatches = Vector{LongDNASeq}()
    for (idx, base) in enumerate(gRNA)
        for b in [DNA_A, DNA_C, DNA_G, DNA_T]
            if b != base
                newgRNA = copy(gRNA)
                newgRNA[idx] = b
                for (idx2, base2) in enumerate(newgRNA)
                    if idx2 == idx
                        continue
                    else
                        for b2 in [DNA_A, DNA_C, DNA_G, DNA_T]
                            if b2 != base2
                                newgRNA2 = copy(newgRNA)
                                newgRNA2[idx2] = b2
                                if newgRNA2 in TwoMismatches
                                    continue
                                else
                                    push!(TwoMismatches, newgRNA2)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    TwoMismatches
end
function readCoordinates(filename::String)
    coordinates = Dict{String,Vector{NamedTuple{(:name, :gRNA, :start, :stop, :strand),Tuple{String,String,Int64,Int64,String}}}}()
    for line in eachline(filename)
        name, gRNA, chrom, start, stop, strand = split(line, '\t')
        if haskey(coordinates, chrom) == false
            coordinates[chrom] = Vector{NamedTuple{(:name, :gRNA, :start, :stop, :strand),Tuple{String,String,Int64,Int64,String}}}()
        end
        push!(coordinates[chrom], (name=name, gRNA=gRNA, start=parse(Int64, start), stop=parse(Int64, stop), strand=strand))
    end
    coordinates
end
function getTail(coord::Dict{String,Vector{NamedTuple{(:name, :gRNA, :start, :stop, :strand),Tuple{String,String,Int64,Int64,String}}}}, fastafile::String, rightFlank::Int64)
    readers = open(FASTA.Reader, fastafile, index=string(fastafile, ".fai"))
    seqs = Dict{String,String}()
    start = 0
    stop = 0
    expand = dna"NNN"
    for (chrom, regions) in coord
        totalChromSequence = FASTA.sequence(readers[chrom])
        for region in regions
            if region.strand == "+"
                start = region.stop + 1
                stop = region.stop + rightFlank
                expand = totalChromSequence[start:stop]
            else
                start = region.start - rightFlank
                stop = region.start - 1
                expand = reverse_complement(totalChromSequence[start:stop])
            end
            if start < 1 || stop > FASTA.seqlen(readers[chrom])
                println("skip $chrom:$(region.start)-$(region.stop) beacuse it is out of range!")
            else
                seqs[region.name] = convert(String, expand)
            end
        end
    end
    seqs
end

function makeBarcode(length::Int64, size::Int64)
    barcodes = BioSequence[]
    i = 1
    while i <= size
        tmp = randdnaseq(length)
        if tmp in barcodes
            continue
        else
            push!(barcodes, tmp)
            i += 1
        end
    end
    barcodes
end

function predictTotalNumbers(gRNA_numbers::Int64)
    OneMis = 3 * 20
    TwoMis = (20 - 1) * 3 * OneMis / 2
    total = OneMis + TwoMis
    Int64(total * gRNA_numbers)
end

@main function makeArtLibrary(; input::String="none", output::String="none", reference::String="none", rightFlank=4, recog="ACCACGTCTCACACC", scaffold="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT", bind="GTTTGGAGACGACGG")
    f = open(output, "w")
    coord = readCoordinates(input)
    tails = getTail(coord, reference, rightFlank)
    barcodes = makeBarcode(8, predictTotalNumbers(length(tails)))
    i = 1
    for (chrom, regions) in coord
        for region in regions
            spacer = LongDNASeq(region.gRNA[1:20])
            totalMismatch = vcat(makeOneMismatch(spacer), makeTwoMismatch(spacer))
            if haskey(tails, region.name)
                tl = tails[region.name]
                for mismatch in totalMismatch
                    chipSequence = string(recog, "g", region.gRNA[1:20], scaffold, "AC", barcodes[i], mismatch, region.gRNA[21:23], tl, bind)
                    println(f, region.name, "\t", chipSequence)
                    i += 1
                end
            end
        end
    end
    close(f)
end
