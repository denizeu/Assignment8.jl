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
        @test uniqueKmers("ACT", 2) == Set(["AC", "CT"])
        @test uniqueKmers("ATC", 2) == Set(["AT", "TC"])
        @test uniqueKmers("ATGCGATG", 4) ==  Set(["TGCG", "ATGC", "GATG", "CGAT", "GCGA"])
        end

    @testset "kmerdist" begin
        @test kmerdist(uniqueKmers("GCGCAT",2), uniqueKmers("ATAT",2)) == 0.8
        @test kmerdist(uniqueKmers("ATCGATG",2), uniqueKmers("GCATACC",2)) == 0.9
    end

    @testset "kmertime" begin
        testpath = normpath(joinpath(@__DIR__, "data"))
        genomes = joinpath(testpath, "cov2_genomes.fasta")
        ex1_path = joinpath(testpath, "ex1.fasta")
        ex2_path = joinpath(testpath, "ex2.fasta")
    

        ex1 = kmertime(parse_fasta(ex1_path)[1], parse_fasta(ex1_path)[2], 3)
        @test ex1 isa Tuple
        @test all(x-> x isa String, ex1[1])
        @test all(x-> x isa String, ex1[2])


    end

    @testset "kmertimes" begin
        path= "AT", "GC", "AC"
        @test kmertimes(path) == 3
    end

    @testset "kmerloc" begin
        Tu= ["ACTA", "AGCT"]
        Japan= ["AGAT", "CCTA", "CCTG"]
        @test kmerloc(path)= Tu=2, Japan=3
    end
end
end

