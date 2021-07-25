using BioinformaticsBISC195
using Test

@testset "BioinformaticsBISC195" begin
    
@testset "Using Strings" begin
    
    @testset "normalizeDNA" begin
        @test normalizeDNA("aatgn") == "AATGN"
        @test_throws Exception normalizeDNA("ZCA")
        @test_throws Exception normalizeDNA(42)
        c = normalizeDNA('C') 
        @test c == "C"
        @test typeof(c) == String
    end # normalizeDNA

    @testset "composition" begin
        seq = rand(['A','T','G','C','N'], 20) |> join
        bc = composition(seq)
        @test bc isa Dict

        @test bc['A'] == count(x-> x == 'A', seq)
        @test bc['C'] == count(x-> x == 'C', seq)
        @test bc['G'] == count(x-> x == 'G', seq)
        @test bc['T'] == count(x-> x == 'T', seq)
        @test bc['N'] == count(x-> x == 'N', seq)

        bc = composition(lowercase(seq))

        @test bc['A'] == count(x-> x == 'A', seq)
        @test bc['C'] == count(x-> x == 'C', seq)
        @test bc['G'] == count(x-> x == 'G', seq)
        @test bc['T'] == count(x-> x == 'T', seq)
        @test bc['N'] == count(x-> x == 'N', seq)
    end # composition

    @testset "gc_content" begin
        @test gc_content("ANTGA") == 0.25 #update method to only count valid bases 
        @test gc_content("cccggg") * 100 == 100.0
        @test gc_content("ATta") == 0.0
        @test_throws Exception gc_content("ATtU")
    end # gc_content

    @testset "complement" begin
        @test complement("ATTAN") == "TAATN"
        @test complement("gcta") == "CGAT"
        @test complement("nnnnnnn") == "NNNNNNN"
        @test_throws Exception complement("APP")
    end # complement

    @testset "reverse_complement" begin
        @test reverse_complement("ATTAN") == "NTAAT"
        @test reverse_complement("gcta") == "TAGC"
        @test reverse_complement("nnnnnnn") == "NNNNNNN"
        @test_throws Exception reverse_complement("AEC")
    end # reverse_complement

    @testset "parse_fasta" begin
        testpath = normpath(joinpath(@__DIR__, "data"))
        genomes = joinpath(testpath, "cov2_genomes.fasta")
        ex1_path = joinpath(testpath, "ex1.fasta")
        ex2_path = joinpath(testpath, "ex2.fasta")

        ex1 = parse_fasta(ex1_path)
        @test ex1 isa Tuple
        @test all(x-> x isa String, ex1[1])
        @test all(x-> x isa String, ex1[2])

        @test ex1[1] == ["ex1.1 | easy", "ex1.2 | multiline"]
        @test ex1[2] == ["AATTATAGC", "CGCCCCCCAGTCGGATT"]

        @test_throws Exception parse_fasta(ex2_path)

        cov2 = parse_fasta(genomes)
        @test length(cov2[1]) == 8
        @test length(cov2[2]) == 8
    end #parse_fasta

    @testset "uniqueKmers" begin
        @test uniqueKmers("GCGCAT",2) isa Set{String}

        @test uniqueKmers("ACT", 2) == Set(["AC", "CT"])
        @test uniqueKmers("ATC", 2) == Set(["AT", "TC"])
        @test uniqueKmers("ATGCGATG", 4) ==  Set(["TGCG", "ATGC", "GATG", "CGAT", "GCGA"])

        @test_throws Exception uniqueKmers("XAQ", 2)
    end

    @testset "kmerdist" begin
        @test kmerdist(uniqueKmers("GCGCAT",2), uniqueKmers("ATAT",2)) isa Float64

        @test kmerdist(uniqueKmers("GCGCAT",2), uniqueKmers("ATAT",2)) == 0.8
        @test kmerdist(uniqueKmers("ATCGATG",2), uniqueKmers("GCATACC",2)) == 0.9

        @test_throws Exception kmerdist(uniqueKmers("AZCGLTG",2), uniqueKmers("GCATACC",2))
    end

    @testset "kmertime" begin
        testpath = normpath(joinpath("BioinformaticsBISC195.jl", "data"))
        genomes = joinpath(testpath, "cov2_genomes.fasta")
        ex1_path = joinpath(testpath, "ex1.fasta")
        ex2_path = joinpath(testpath, "ex2.fasta")
        datatry_path= joinpath(testpath, "datatry.fasta")
    

        ex1 = kmertime(parse_fasta(ex1_path)[1], parse_fasta(ex1_path)[2], 3)
        @test ex1 isa Tuple
        @test all(x-> x isa String, ex1[1])
        @test all(x-> x isa String, ex1[2])

        datatry= kmertime(parse_fasta(datatry_path)[1], parse_fasta(datatry_path)[2], 2)
        @test datatry== (Any[Set(["AC", "CT"])], Any[Set(["CC", "GC", "GG", "CG", "AT", "CA", "TG", "TA", "GT", "GA", "TT", "AC", "CT", "AA", "AG", "TC"])], Any[Set(["AG", "GG", "GT"])]) 
    end

    @testset "kmertimes" begin
        testpath = normpath(joinpath(@__DIR__, "..", "data"))
        ex1_path = joinpath(testpath, "ex1.fasta")
        headers, sequences = parse_fasta("data/refined_data.fasta");

        @test kmertimes(headers, sequences) isa Tuple{Int64, Int64, Int64}
        @test kmertimes(headers, sequences) == (64, 73, 95)
        @test_throws Exception kmertimes(parse_fasta(ex1_path)[1], parse_fasta(ex1_path)[2])
    end

    @testset "pairdist" begin
        testpath = normpath(joinpath(@__DIR__, "data"))
        ex2_path = joinpath(testpath, "ex2.fasta")
        refine_path = joinpath(testpath,"refined_data.fasta")

        @test pairdist(refine_path) isa Matrix{Float64}
        @test pairdist(refine_path)[6] == 0.0072727272727273196
        @test_throws Exception pairdist(ex2_path)
        
    end

    @testset "distsort" begin
        testpath = normpath(joinpath(@__DIR__, "data"))
        refine_path = pairdist(joinpath(testpath, "refined_data.fasta"))
        datasort_path = joinpath(testpath, "datasort.fasta")

        
        @test distsort(refine_path) isa Tuple{Vector{Any}, Vector{Any}, Vector{Any}}
        @test distsort(refine_path)[1][121] == 0.006805257760790551
        @test_throws Exception distsort(datasort_path)
    end
end
end

