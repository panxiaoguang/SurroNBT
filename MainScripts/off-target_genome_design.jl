using Fire
using FASTX
using BioSequences
## make a struct to store the off-target information
struct OffTarget
    ont::AbstractString
    chrom::AbstractString
    startSite::Int64
    offSeq::AbstractString
    strand::AbstractString
    mismatches::Int64
end
## make input
function makeInput(targetSeq::Vector{String}; ref_dir::String="/home/panxiaoguang/Project/DataBase/GRCh38_chrom", spacer_length::Int64=20, PAM::String="NRG", MaxMismatches::Int64=4)
    (path, io) = mktemp()
    println(io, ref_dir)
    println(io, repeat('N', spacer_length), PAM)
    for tg in targetSeq
        println(io, tg, " ", MaxMismatches)
    end
    close(io)
    path
end
## run cas-offinder
function runCasoffinder(; targetSeq::String="none", output::String="none", ref_dir::String="/home/panxiaoguang/Project/DataBase/GRCh38_chrom", spacer_length::Int64=20, PAM::String="NRG", MaxMismatches::Int64=4)
    targets = Vector{String}()
    for line in eachline(targetSeq)
        push!(targets, line)
    end
    inputs = makeInput(targets; ref_dir=ref_dir, spacer_length=spacer_length, PAM=PAM, MaxMismatches=MaxMismatches)
    run(`cas-offinder $inputs C $output`)
end

## parse output
function parseOutput(output::String)
    offTargets = Vector{OffTarget}()
    for line in eachline(output)
        ont, chrom, startSite, offSeq, strand, mismatches = split(line, '\t')
        startSite = parse(Int64, startSite)
        mismatches = parse(Int64, mismatches)
        push!(offTargets, OffTarget(ont, chrom, startSite, offSeq, strand, mismatches))
    end
    offTargets
end

## main function

function getOffs(; targetSeq::String="none", ref_dir::String="/home/panxiaoguang/Project/DataBase/GRCh38_chrom", spacer_length::Int64=20, PAM::String="NRG", MaxMismatches::Int64=4)
    result = tempname()
    runCasoffinder(targetSeq=targetSeq, output=result, ref_dir=ref_dir, spacer_length=spacer_length, PAM=PAM, MaxMismatches=MaxMismatches)
    finallyResult = parseOutput(result)
    finallyResult
end

function getTail(coord::Vector{OffTarget}, fastafile::String, rightFlank::Int64)
    readers = open(FASTA.Reader, fastafile, index=string(fastafile, ".fai"))
    seqs = Dict{String,String}()
    start = 0
    stop = 0
    expand = dna"NNN"
    for off in coord
        if off.strand == "+"
            start = off.startSite + 1
            stop = off.startSite + 23 + rightFlank
            expand = FASTA.extract(readers, DNAAlphabet{4}(), off.chrom, start:stop)
        else
            start = off.startSite + 1 - rightFlank
            stop = off.startSite + 23
            expand = reverse_complement(FASTA.extract(readers, DNAAlphabet{4}(), off.chrom, start:stop))
        end
        if start < 1 || stop > FASTA.seqlen(readers[off.chrom])
            println("skip $(off.chrom):$(off.startSite) beacuse it is out of range!")
        else
            seqs["$(off.chrom):$(off.startSite)"] = convert(String, expand)
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


@main function makeGenLibrary(; input::String="none", output::String="none", ref_dir::String="/home/panxiaoguang/Project/DataBase/GRCh38_chrom", spacer_length::Int64=20, PAM::String="NRG", MaxMismatches::Int64=4, reference::String="none", rightFlank=4, recog="ACCACGTCTCACACC", scaffold="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT", bind="GTTTGGAGACGACGG")
    f = open(output, "w")
    offs = getOffs(targetSeq=input, ref_dir=ref_dir, spacer_length=spacer_length, PAM=PAM, MaxMismatches=MaxMismatches)
    tails = getTail(offs, reference, rightFlank)
    barcodes = makeBarcode(8, length(tails))
    i = 1
    for off in offs
        gRNA = off.ont[1:20]
        name = "$(off.chrom):$(off.startSite)"
        if haskey(tails, name)
            tl = tails[name]
            chipSequence = string(recog, "g", gRNA, scaffold, "AC", barcodes[i], tl, bind)
            println(f, off.chrom, "\t", off.startSite, "\t", off.startSite + 23, "\t", off.offSeq, "\t", off.strand, "\t", off.mismatches, "\t", chipSequence)
            i += 1
        end
    end
    close(f)
end

