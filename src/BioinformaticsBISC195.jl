module BioinformaticsBISC195

using Base: stop_reading
export normalizeDNA
export composition
export gc_content
export complement
export reverse_complement
export parse_fasta
export uniqueKmers
export kmerdist
export kmerloc
export kmertime
export kmertimes
export lengthcount
export minMax
export pairdist
export distsort

# ### 1. NormalizeDNA Function
"""
    normalizeDNA(::AbstractString)

Function that takes a sequence and returns the uppercase and normalized version of it.

Returns an error if invalid bases are encountered.

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> normalizeDNA("acatg")
    "ACATG"

    julia> normalizeDNA("LAC")
    ERROR: key 'L' not found

    julia> normalizeDNA(7)
    ERROR: MethodError: no method matching uppercase(::Int64)

    julia> normalizeDNA('G')
    "G"

    julia> s = normalizeDNA('C');

    julia> typeof(s)
    String
"""
function normalizeDNA(sequence)
    sequence= uppercase(sequence)
    rep = Dict('R' => 'N', 'Y' => 'N', 'S' => 'N', 'K' => 'N', 'M' => 'N', 'B' => 'N', 'D' => 'N', 'H' => 'N', 'W' => 'N','V' => 'N', '.' => 'N', '-' => 'N', 'A' => 'A', 'G' => 'G', 'C' => 'C', 'T' => 'T', 'N' => 'N')
    return join([rep[c] for c in sequence])
end

# ### 2. Composition Function
"""
    composition(sequence)

Takes a DNA sequence and counts the number of each type of base in the DNA sequence.

Returns a Dictionary of how many of each exist in the order A, C, G, T.

Converts any ambiguous bases not in AGCT to "N".

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> composition("ACGGG")
    'A' => 1
    'G' => 3
    'C' => 1

    julia> composition("ACC")
    'A' => 1
    'C' => 2

    julia> A,C,T = composition("acctt")
    'A' => 1
    'T' => 2
    'C' => 2

    julia> composition("BAGGGRRRR")
    'A' => 1
    'G' => 3
    'N' => 5   
"""
function composition(sequence)
    sequence = normalizeDNA(sequence) #make uppercase string, check invalid bases
    comp = Dict()
    for i in 1:length(sequence)
        pos = sequence[i]
        if haskey(comp, pos) == true
            comp[pos] = comp[pos] + 1
        else
            comp[pos] = 1
        end
    end
    return comp
end

# ### 3. GC_Content Function
"""
    gc_content(sequence)

Calculates the GC ratio of a DNA sequence.

The GC ratio is the total number of G and C bases divided by the total number of G's, C's, A's, and T's within the sequence.

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> gc_content("CATA")
    0.25

    julia> gc_content("cccg") * 20
    20.0

    julia> gc_content("AaTt")
    0.0

    julia> gc_content("ATty")
    0.0
    
    julia> gc_content("GGCCNNNN")
    1.0
"""
function gc_content(sequence)
    comp = composition(sequence)
    ng = get(comp, 'G', 0) #searches for each of these bases using the composition() function above
    nc = get(comp, 'C', 0)
    na= get(comp, 'A', 0)
    nt= get(comp, 'T', 0)
    return (ng + nc) / (ng + nc + nt + na) #divided by sum of AGCT w/ composition function
end

# ###4. Complement Function
"""
    complement(base)

Get the DNA complement of the provided base:

    A <-> T
    G <-> C

Accepts uppercase or lowercase `String` or `Char`,
but always returns an uppercase `String` (orig says Char).

If a valid base is not provided, the function returns an error.

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> complement("A")
    "T"

    julia> complement("GCT")
    "CGA"

    julia> complement("AGCTACC")
    "TCGATGG"

    julia> complement("ZAG")
    ERROR: KeyError: key 'Z' not found

"""
function complement(sequence)
    sequence = normalizeDNA(sequence)
    ret = "" #initialize empty string to place the complements
    for i in 1:length(sequence)
        if sequence[i] == 'A'
            ret = ret * 'T' #if position equals A return "T" in that new position as the complement
        elseif sequence[i] == 'T'
            ret = ret * 'A'
        elseif sequence[i] == 'G'
            ret = ret * 'C'
        elseif sequence[i] == 'C'
            ret = ret * 'G'
        elseif sequence[i] == 'N'
            ret = ret * 'N'
        end
    end
    return ret 
end

# ###5. Reverse Complement Function
"""
    reverse_complement(sequence)

Takes a DNA sequence and returns the reverse complement
of that sequence.
    
Takes lowercase or uppercase sequences, but always returns uppercase.
    
    Examples
    ≡≡≡≡≡≡≡≡≡≡
        julia> reverse_complement("ACCTTT")
        "AAAGGT"
    
        julia> reverse_complement("CCGTAGTA")
        "TACTACGG"

        julia> reverse_complement("ttacg")
        "CGTAA"
    
        julia> rc = reverse_complement("ACAG");
    
        julia> println(rc)
        CTGT
    """
function reverse_complement(sequence)
    st= "" #initialize empty string to place the reverse complements
    sequence =reverse(sequence) 
    for i in 1:length(sequence)
        st = st*string(complement(sequence[i])) 
    end
    return st
end

# ### Parse Function

"""
    function parse_fasta(path)

Reads a fasta-formated file and returns 2 vectors,
one containing parsed headers,
the other containing the entire sequence as a `String`.

Note: function does not validate DNA sequences for correctness.

Example
≡≡≡≡≡≡≡≡≡
    julia> ex1 = parse_fasta("data/ex1.fasta");

    julia> ex1 isa Tuple
    true

    julia> ex1[1]
    2-element Array{Tuple{String,String},1}:
      ("ex1.1", "easy")
      ("ex1.2", "multiline")

    julia> ex1[2]
    2-element Array{String,1}:
    "AATTATAGC"
    "CGCCCCCCAGTCGGATT"

    julia> ex2 = parse_fasta("data/ex2.fasta");

    julia> ex2[1]
    4-element Array{Tuple{String,String},1}:
      ("ex2.1", "oneper")
      ("ex2.2", "wrong")
      ("ex2.3", "wronger")
      ("ex2.4", "wrongest")
    
    julia> ex2[2]
    4-element Array{String,1}:
      "ATCCGT"
      "ATCGTGGaact"
      "ATCGTGGaact"
      "this isn't a dna string,but parse it anyway"
"""
function parse_fasta(path)
    headers= String[] #initialize empty array to hold the headers
    sequences= [] #initialize empty array to hold sequences
    tmp= "" 
    for line in eachline(path)
        if startswith(line, '>')
            if !isempty(tmp)
                push!(sequences, tmp) #if it is empty, put the sequences within tmp
            end
            tmp= ""
            push!(headers, SubString(line, 2))
        else
            line = normalizeDNA(line)
            tmp= tmp*line
        end
    end
        push!(sequences, tmp)
    return (headers, sequences)
end

# ### Mean Lengths and Counts
"""
    function lengthcount(path)

Takes data and returns the mean and standard deviation of sequence lengths and of gc content within one tuple.
    
    
    Example
    ≡≡≡≡≡≡≡≡≡≡
        julia> using Statistics
    
        julia> lengthcount("data/datatry.fasta")
        (23.666666666666668, 34.93326972004386, 0.3975694444444444, 0.08965799537274553)

        julia> lengthcount("data/ex1.fasta")
        (13.0, 5.656854249492381, 0.46405228758169936, 0.34199935822094457)
    """
function lengthcount(path)
    lengths= []
    counts= []
    data= parse_fasta(path) #parsing the data in order to separate headers and sequences
    for i in data[2] 
        push!(lengths, length(i)) 
        gs = count(==('G'), i) #counts the number of g's
        cs = count(==('C'), i) #counts the number of c's
        GCcontent= (gs + cs)/length(i) #finds the GC content by dividing by length of sequence
        push!(counts, GCcontent)
    end
    return (mean(lengths), std(lengths), mean(counts), std(counts)) #using Statistics, takes mean and SD of each calculation
end

# ### Minimum and Maximum Function
"""
    function minMax(path)

Takes data and returns a Tuple of the minimum and maximum of sequence lengths.
    
    Examples
    ≡≡≡≡≡≡≡≡≡≡
        julia> minMax("data/ex1.fasta")
        (9, 17)

        julia> minMax("data/datatry.fasta")
        (3, 64)

        julia> minMax("data/datasort.fasta")
        (21, 51)

        julia> typeof(minMax("data/datasort.fasta"))
        Tuple{Int64, Int64}
    """
function minMax(path)
    ret = parse_fasta(path)[2]
    retmin = length(minimum(ret))
    retmax = length(maximum(ret))
    return (retmin, retmax)
end

# ### Sequence Lengths Function
"""
        seqlength(path)

Takes data and returns a vector of the sequence lengths within the data.

    Examples
    ≡≡≡≡≡≡≡≡≡≡
        julia> seqlength("data/ex1.fasta")
        2-element Vector{Any}:
            9
            17

        julia> seqlength("data/datatry.fasta")
        3-element Vector{Any}:
            3
            64
            4

    """
function seqlength(path)
    lengths = []
    data = parse_fasta(path)
    for i in data[2]
         push!(lengths, length(i))
    end
    return lengths
end

# ### Unique Kmer Function
"""
    uniqueKmers(sequence, k)

Takes a sequence and a kmer length (k) 
and returns a set of strings of the unique kmers that appear within the DNA sequence.

Converts invalid bases to "N" and converts lowercase sequences to uppercase sequences.

Example
≡≡≡≡≡≡≡≡≡
    julia> uniqueKmers("ATGCGATG", 4)
    Set{Any} with 5 elements:
    "TGCG"
    "ATGC"
    "GATG"
    "CGAT"
    "GCGA"

    julia> uniqueKmers("ATGCN", 2)
    Set{String} with 4 elements:
    "CN"
    "AT"
    "GC"
    "TG"

    julia> uniqueKmers("ACT", 2)
    Set{Any} with 2 elements:
    "AC"
    "CT"

    uniqueKmers("RAcGcN", 4)
    Set{String} with 3 elements:
    "NACG"
    "CGCN"
    "ACGC"
"""
function uniqueKmers(sequence, k)
    sequence= normalizeDNA(sequence)
    1 <= k <= length(sequence) || error("k must be a positive integer less than the length of the sequence")
    kmers = Set{String}()    
    stopindex = length(sequence) - k + 1
    for i in 1:stopindex
        kmer= sequence[i:i+k-1]
        kmer = uppercase(kmer)
        push!(kmers, kmer)
    end
    return kmers
end

# ### Sorting Sequences to 5%
"""
    function sortingseq(path)

Takes a path (data set) and returns vectors:
one with the sequence lengths greater than 29500 and one with the corresponding headers.

Example
≡≡≡≡≡≡≡≡≡
    julia> sorting("data/data_sorting.fasta")
    (Any[29903, 29838, 29771, 29782], ["NC_045512.2 |China|Homo sapiens", "MT614596.1 |Egypt|Homo sapiens", "MT669321.1 |India|Homo sapiens", "MT670022.1 |Chile|Homo sapiens"])

    julia> g= sorting("data/data_sorting.fasta");

    julia> typeof(g)
    Vector{Any} (alias for Array{Any, 1})

"""
function sortingseq(path)
    sequences = seqlength(path)
    headers = parse_fasta(path)[1]
    check = findall(x->x<29500, sequences)
    deleteat!(sequences, check)
    deleteat!(headers, check)
    return (sequences, headers)
end

##Histogram for lengths: see analysis repo
#=data= sorting("data/refined_data.fasta")[1];
histogram(data, bins= 10, label= "Sorted Sequences", xlabel= "Sequence Lengths", ylabel= "Number of Sequences", legend=:topleft)=#

# ### 8. KmerDistance Function
"""
    function kmerdist(set1, set2)

Takes two kmer sets and returns the distance between the kmer sets as a Float64.
Must return a positive number between 0 and 1.

Returns 0.0 for identical sets.
Returns 1.0 for sets with no similarity.

Example
≡≡≡≡≡≡≡≡≡
    julia> kmerdist(uniqueKmers("ATCGATG",2), uniqueKmers("GCATACC",2))
    0.9

    julia> kmerdist(uniqueKmers("GCGCAT",2), uniqueKmers("ATAT",2))
    0.8

    julia> kmerdist(uniqueKmers("AATA", 4), uniqueKmers("CGGCCCG", 4))
    1.0

    julia> typeof(kmerdist(uniqueKmers("ATCGATG",2), uniqueKmers("GCATACC",2)))
    Float64
"""
function kmerdist(set1, set2)
    return 1 - (length(intersect(set1, set2))/length(union(set1,set2))) #length of intersection of the sets divided by length of the union
end

# ###9. Kmertime Function
"""
    function kmertime(headers, sequences, k)

Takes a dataset and returns three vectors:
of the unique kmers along with corresponding headers for early, middle, and late time periods.

In this data: (and my project: "data/refined_data.fasta"), early: 2019, middle: 2020, late: 2021.

Example
≡≡≡≡≡≡≡≡≡
    julia> headers, sequences = parse_fasta("data/datatry.fasta");

    julia> kmertime(headers, sequences, 3)
    (Any[Set(["ACT"])], Any[Set(["TGT", "CTT", "GAC", "AGT", "GAA", "TGC", "TTC", "ACC", "ACA", "TAG"  …  "CCA", "TTT", "AGA", "TGG", "ATA", "AAG", "CAC", "AAA", "TTA", "CCG"])], Any[Set(["GGT", "AGG"])])

    julia> j= kmertime(headers, sequences, 3);

    julia> typeof(j)
    Tuple{Vector{Any}, Vector{Any}, Vector{Any}}
"""
function kmertime(headers, sequences, k=3) #whatto set k to?
    early= []
    middle= []
    late= [] #Initializing empty arrays for 3 time periods: early(2019), middle(2020), and late(2021)
    for i in 1:length(sequences)
        if occursin("2019", headers[i]) #Findall occurrences of 2019 per header
            push!(early, uniqueKmers(sequences[i], k)) #If the header contains the date "2019", the kmer will be pushed into the "early" array, calling the uniqueKmer function to process how many unique kmers exist in the sequence
        end
        if occursin("2020", headers[i])
            push!(middle, uniqueKmers(sequences[i], k)) #If the header contains the date "2020", the kmer will be pushed into the "middle" array
        end
        if occursin("2021", headers[i])
            push!(late, uniqueKmers(sequences[i], k)) #If the header contains the date "2021", the kmer will be pushed into the "late" array
        end
    end
    return (early, middle, late)
end


# ###Kmertimes Function
"""
    function kmertimes(path)

Takes a dataset and returns an Tuple of integers for 
the number of unique kmers within early, middle, and late time periods.

Example
≡≡≡≡≡≡≡≡≡

    julia> headers, sequences = parse_fasta("data/datatry.fasta");

    julia> kmertimes(headers, sequences)

    julia> (1, 36, 2)

    julia> typeof(kmertimes(headers, sequences))
    Tuple{Int64, Int64, Int64}
"""
function kmertimes(headers, sequences)
    early, middle, late = kmertime(headers, sequences) #calling kmertime to separate the sequences by time period: 2019, 2020, 2021
    early_kmers = length(union(early...)) #calculating length of union of the sets within each time period
    middle_kmers = length(union(middle...))
    late_kmers = length(union(late...))
    return(early_kmers, middle_kmers, late_kmers) #returns the number of unique kmers in each time period
end

### Kmertimes Plot: See Analysis Repo
"""
    Returns a barplot for kmertime: returns one plot for the number of unique kmers in early middle and late time periods. 

    The time periods are represented by bars of 3 different colors corresponding to each different time period.
"""
#=bar(["early" "middle" "late"],
           [64 73 95],
           labels = ["early" "middle" "late"],
           label = "Number of Unique Kmers",
           title = "Time Period vs. Number of Unique Kmers",
           xlabel = "Time Period",
           ylabel = "Number of Unique Kmers",
           color = [:steelblue :pink :lavender],
           bg= "beige",
           legend = :topleft)=#

# ###10. Pairwise Distance Function
"""
    pairdist(path)
    
Takes a dataset and initiates a matrix for sequences and
returns the distance between pairs of sequences.

Example
≡≡≡≡≡≡≡≡≡
    julia> ("data/refined_data.fasta")[35][1]
    0.189873417721519

    julia> pairdist("data/refined_data.fasta")[70][1]
    0.24705882352941178

    julia> typeof(pairdist("data/refined_data.fasta"))
    Matrix{Float64} (alias for Array{Float64, 2})
"""

function pairdist(path) 
    h = kmertime(parse_fasta(path)[1], parse_fasta(path)[2], 3)[1] #takes kmertime of headers and sequences and stores them with these letters
    j =  kmertime(parse_fasta(path)[1], parse_fasta(path)[2], 3)[2]
    k =  kmertime(parse_fasta(path)[1], parse_fasta(path)[2], 3)[3]
    mesh = vcat(h,j,k) #concatenates h, j, k
    ret = zeros(36, 36) #initiates empty matrix
    for i in 1:36 
        for j in 1:36
            i <= j && continue #skips upper triangle, focuses bottom half of matrix
            d = kmerdist(mesh[i], mesh[j]) #finds distance between the kmers
            ret[i, j] = d
        end
    end
    return ret #returns the matrix with the added distance between kmers of sequences
end

### Distance Sort
"""
    distsort(mat)
    
Takes a matrix (of my "data/refined_data.fasta") and returns a Tuple of Vectors 
containing the distances between unique kmers of the sequences, organized by time
period.

Example
≡≡≡≡≡≡≡≡≡
    julia> distsort(pairdist("data/refined_data.fasta"))[3][142]
    0.16883116883116878

    julia> typeof(distsort(pairdist("data/refined_data.fasta")))
    Tuple{Vector{Any}, Vector{Any}, Vector{Any}}
"""

function distsort(mat) #takes matrix as an argument
    early_dist = [] #3 separate one dimensional arrays initialized for each time period
    middle_dist = []
    late_dist = []
    for i in 1:12
        for j in 1:12 #goes through the early time periods (1:12 are the first 12 in my refined_data)
            push!(early_dist, mat[i,j]) 
        end
    end
    for i in 13:24 #goes thru the middle time periods
        for j in 13:24
            push!(middle_dist, mat[i,j])
        end
    end
    for i in 25:36 #goes thru the late time periods
        for j in 25:36
            push!(late_dist, mat[i,j])
        end
    end
    return (early_dist, middle_dist, late_dist)
end
end

        

### Sorted Distance BoxPlot
#=early_dist= distsort(pairdist("data/refined_data.fasta"))[1]
middle_dist =distsort(pairdist("data/refined_data.fasta"))[2]
late_dist = distsort(pairdist("data/refined_data.fasta"))[3]
boxplot(["early", "middle", "late"], [[early_dist], [middle_dist], [late_dist]]
       , title="Distances vs Time Periods BoxPlot", xlabel="Distances", ylabel="Time Periods")  =#

