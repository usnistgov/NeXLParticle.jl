using CSV
using DataStructures: SortedDict, SortedSet
using Random
using DataAPI
using StringEncodings
using StaticArrays: SVector
using CoordinateTransformations: AbstractAffineMap

Base.convert(::Type{Symbol}, elm::Element) = Symbol(uppercase(elm.symbol))


"""
Handles RJ Lee-style Zeppelin particle data sets.
"""
struct Zeppelin
    headerfile::String
    header::SortedDict{String,String}
    data::DataFrame
    classnames::Vector{String}
    
    function Zeppelin( #
        headerfile::String,
        header::AbstractDict{String,String},
        data::DataFrame,
        classnames::AbstractVector{String}
    )
        function _massagehdz(header)
            res = SortedDict(header)
            for hi in keys(header)
                if  hi in ( "PARTICLE_PARAMETERS", "PARAMETERS", "HEADER_FMT" ) || # 
                    startswith(hi, "ELEM") || startswith(hi, "CLASS")
                    delete!(res, hi)
                end
            end
            return res
        end
        # Mostly converts to pretty printing
        function _massagepxz(zep)
            zd = zep.data
            # Find all the class name mappings
            for col in ("CLASS", "VERIFIEDCLASS")
                ic = findfirst(nm->nm == col, names(zd))
                if (!isnothing(ic)) && (eltype(zd[:,ic]) <: Integer)
                    old = zd[:, ic]
                    select!(zd, Not(ic))
                    insertcols!(zd, ic, col => map(ncl->ZepClass(zep, ncl), old))
                end
            end
            for col in ( "FIRSTELM", "SECONDELM", "THIRDELM", "FOURTHELM" )
                ic = findfirst(nm->nm==col, names(zd))
                if !isnothing(ic) && (eltype(zd[:,ic]) <: Integer)
                    old = zd[:, ic]
                    select!(zd, Not(ic))
                    insertcols!(zd, ic, col=>map(z->z in 1:100 ? elements[z] : missing, old))
                end
            end
            for col in ( "XDAC", "YDAC", "TYPE_4ET_", "TYPE4ET" )
                ic = findfirst(nm->nm==col, names(zd))
                if !isnothing(ic)
                    select!(zd, Not(ic))
                end
            end
        end
        zep = new(
            endswith(headerfile, r"\.[h|H][d|D][z|Z]$") ? headerfile : "$headerfile.hdz",
            _massagehdz(header),
            copy(data),
            classnames
        )
        _massagepxz(zep)
        return zep
    end

    Zeppelin( hdzfilename::String) = loadZep( hdzfilename)
end


function NeXLSpectrum.name(zep::Zeppelin)
    get(zep.header, "NAME", get(zep.header, "DESCRIPTION", splitpath(zep.headerfile)[end-1]))
end


struct ZepClass
    zep::Zeppelin
    index::Int
end
function Base.show(io::IO, zc::ZepClass)
    if zc.index == -1
        print(io, "-")
    elseif zc.index+1 in eachindex(zc.zep.classnames)
        print(io, zc.zep.classnames[zc.index+1])
    else
        print(io, "**UNCLASSIFIED**")
    end
end
Base.:(==)(zc1::ZepClass, zc2::ZepClass) = zc1.index==zc2.index
Base.:(==)(zc1::ZepClass, zc2::AbstractString) = repr(zc1)==zc2
Base.:(==)(zc1::AbstractString, zc2::ZepClass) = zc1==repr(zc2)
Base.String(zc::ZepClass) = repr(zc)

Base.get(z::Zeppelin, key::String, def=missing) = get(z.header, uppercase(key), def)
Base.getindex(z::Zeppelin, key::String) = getindex(z.header, uppercase(key))

Base.copy(z::Zeppelin) = Zeppelin(z.headerfile, copy(z.header), copy(z.data), copy(z.classnames))

function loadZep( hdzfilename::String)::Zeppelin
    remapcolumnnames = Dict(
        "PART#" => "NUMBER",
        "PARTNUM" => "NUMBER",
        "FIELD#" => "FIELD",
        "FIELDNUM" => "FIELD",
        "MAGFIELD#" => "MAGFIELD",
        "MAGFIELDNUM" => "MAGFIELD",
        "X_ABS" => "XABS",
        "Y_ABS" => "YABS",
        "X_DAC" => "XDAC",
        "Y_DAC" => "YDAC",
        "XCENT" => "XDAC",
        "YCENT" => "YDAC",
        "X_FERET" => "XFERET",
        "Y_FERET" => "YFERET",
        "DVG" => "DAVG",
        "DAVE" => "DAVG",
        "PERIM" => "PERIMETER",
        "ORIENT" => "ORIENTATION",
        "LIVE_TIME" => "LIVETIME",
        "FIT_QUAL" => "FITQUAL",
        "MAG_INDEX" => "MAGINDEX",
        "FIRST_ELEM" => "FIRSTELM",
        "SECOND_ELEM" => "SECONDELM",
        "THIRD_ELEM" => "THIRDELM",
        "FOURTH_ELEM" => "FOURTHELM",
        "ATOMICNUMBER1" => "FIRSTELM",
        "ATOMICNUMBER2" => "SECONDELM",
        "ATOMICNUMBER3" => "THIRDELM",
        "ATOMICNUMBER4" => "FOURTHELM",
        "FIRST_CONC" => "COUNTS1",
        "SECOND_CONC" => "COUNTS2",
        "THIRD_CONC" => "COUNTS3",
        "FOURTH_CONC" => "COUNTS4",
        "FIRST_PCT" => "FIRSTPCT",
        "SECOND_PCT" => "SECONDPCT",
        "THIRD_PCT" => "THIRDPCT",
        "FOURTH_PCT" => "FOURTHPCT",
        "PCT1" => "FIRSTPCT",
        "PCT2" => "SECONDPCT",
        "PCT3" => "THIRDPCT",
        "PCT4" => "FOURTHPCT",
        "TYPE(4ET)#" => "TYPE4ET",
        "TYPE(4ET)" => "TYPE4ET",
        "VOID_AREA" => "VOIDAREA",
        "RMS_VIDEO" => "RMSVIDEO",
        "FIT_QUAL" => "FITQUAL",
        "VERIFIED_CLASS" => "VERIFIEDCLASS",
        "EDGE_ROUGHNESS" => "EDGEROUGHNESS",
        "COMP_HASH" => "COMPHASH",
        "PSEM_CLASS" => "CLASS",
        "VERIFIED_CLASS" => "VERIFIEDCLASS",
    )
    function _extractclassnames(header)
        # Find all the class name mappings
        cn = SortedDict(parse(Int, cn[6:end]) => header[cn] for cn in filter(c -> !isnothing(match(r"^CLASS\d+", c)), collect(keys(header))))
        return map(i -> get(cn, i, "CLASS$i"), 0:maximum(keys(cn)))
    end    
    columnnames(cols) = uppercase.(map(cn -> get(remapcolumnnames, cn, cn), map(c -> c[1], cols)))
    header, columns, hdr = Dict{String,String}(), [], true
    open(hdzfilename, enc"WINDOWS-1252","r") do f
        for line in readlines(f)
            if hdr
                p = findfirst(c -> c == '=', line)
                if !isnothing(p)
                    (k, v) = line[1:p-1], line[p+1:end]
                    header[k] = v
                    hdr = !isequal(uppercase(k), "PARTICLE_PARAMETERS")
                end
            else
                push!(columns, string.(strip.(split(line, "\t"))))
            end
        end
    end
    pxz = CSV.File(
        replace(hdzfilename, r".[h|H][d|D][z|Z]$" => ".pxz"),
        header = columnnames(columns),
        delim='\t',
        normalizenames = true,
        missingstring="-",
        decimal='.'
    ) |> DataFrame
    return Zeppelin( hdzfilename, header, pxz, _extractclassnames(header))
end

function Base.show(io::IO, zep::Zeppelin)
    print(io, "Zeppelin[$(zep.headerfile),$(size(zep.data))]")
end

Base.size(z::Zeppelin) = size(z.data)
Base.size(z::Zeppelin, i::Int) = size(z.data, i)

classes(zep::Zeppelin) = levels(zep[:,"CLASS"])

NeXLCore.elms(zep::Zeppelin) = SortedSet(filter(elm -> uppercase(elm.symbol) in names(zep.data), PeriodicTable.elements[1:94]))

header(zep::Zeppelin) = SortedDict(zep.header)

"""
    asa(::Type{DataFrame}, zep::Zeppelin; rows=missing, columns=missing, sortcol=:None, rev=false)

Create a DataFrame containing the particle data.  Optionally specify which rows and columns to include and
whether to sort the resulting table.
"""
function NeXLUncertainties.asa(::Type{DataFrame}, zep::Zeppelin; rows=missing, sortcol=:None, columns=missing, rev=false)
    rs = ismissing(rows) ? eachparticle(zep) : intersect(eachparticle(zep), rows)
    cols = ismissing(columns) ? names(zep.data) : intersect(columns,names(zep.data))
    res = zep.data[rs, cols]
    return sortcol≠:None ? sort(res, sortcol, rev=rev) : res
end

eachparticle(zep::Zeppelin) = 1:size(zep.data, 1)

function DataAPI.describe(zeps::AbstractVector{Zeppelin}; dcol=:DAVG, nelms=2)
    function build(zep)
        df=describe(zep, dcol=dcol, nelms=nelms)
        insertcols!(df, 1, :Dataset=>[zep.header["DESCRIPTION"] for _ in 1:size(df,1)])
        return df
    end
    return vcat(map(z->build(z),zeps))
end

function DataAPI.describe(zep::Zeppelin; dcol=:DAVG, nelms=3)
    els = elms(zep)
    function sortedstats(rows)
        tmp = [ ( elm, summarystats(zep[rows,convert(Symbol, elm)])) for elm in els ]
        return sort!(tmp, lt=(v1,v2)->isless(v1[2].mean,v2[2].mean),rev=true)
    end
    cnelms = min(length(els), nelms)
    ds = summarystats(zep[:,dcol])
    clss = String["All"]
    szs = Int[size(zep.data,1)]
    minds = Float64[ ds.min ]
    medds = Float64[ ds.median ]
    maxds = Float64[ ds.max ]
    sbcm = StatsBase.countmap(zep[:,"CLASS"])
    for (cn, _) in sbcm
        rows = rowsclass(zep,repr(cn))
        #if length(rows)>0
        push!(clss, repr(cn))
        push!(szs, length(rows))
        ds = summarystats(zep[rows,dcol])
        push!(minds, ds.min)
        push!(medds, ds.median)
        push!(maxds, ds.max)
        #end
    end
    dfres = DataFrame(Symbol("Class")=>clss, Symbol("Count")=>szs,Symbol("Min[$dcol]")=>minds, Symbol("Median[$dcol]")=>medds, Symbol("Max[$dcol]")=>maxds)
    # Add All elemental row
    ss = sortedstats(eachparticle(zep))
    elmdfs=DataFrame[ ]
    uem, uef = Union{String,Missing}, Union{Float64,Missing}
    for ne in 1:cnelms
        cols = Symbol.( [ "Elm[$ne]", "Min[$ne]", "Median[$ne]", "Max[$ne]" ])
        sst = ss[ne][2]
        push!(elmdfs, DataFrame(cols[1]=>uem[ symbol(ss[ne][1]) ], cols[2]=>uef[sst.min], cols[3]=>uef[sst.median], cols[4]=>uef[sst.max]))
    end
    for (cn, _) in sbcm
        try
            ss = sortedstats(rowsclass(zep, repr(cn)))
            for ne in 1:cnelms
                sst = ss[ne][2]
                push!(elmdfs[ne], [ symbol(ss[ne][1]), sst.min, sst.median, sst.max ])
            end
        catch
            for ne in 1:cnelms
                push!(elmdfs[ne], [ missing, missing, missing, missing ])
            end
        end
    end
    for elmdf in elmdfs
        for col in names(elmdf)
            insertcols!(dfres, size(dfres,2)+1, col=>elmdf[:,col])
        end
    end
    return sort!(dfres, [:Count, :Class], rev=true)
end


"""
    zep[123] # where zep is a Zeppelin
    zep[1:2:8]

Returns the Spectrum (with images) associated with the particle at row or rows
"""
Base.getindex(zep::Zeppelin, row::Int) = spectrum(zep, row, true)
Base.getindex(zep::Zeppelin, rows) = map(row->spectrum(zep, row, true), rows)
Base.getindex(zep::Zeppelin, rows, cols) = getindex(zep.data, rows, cols)
Base.lastindex(zep::Zeppelin, axis::Integer) = lastindex(zep.data, axis)

# Replace the default in PeriodicTable because it is too verbose...
# Base.show(io::IO, elm::Element) = print(io, elm.symbol)


function filenumber(zep::Zeppelin, row::Int)::String
    tmp="$(zep.data[row,:NUMBER])"
    return repeat('0',max(0,5-length(tmp)))*tmp
end

"""
    spectrumfilename(zep::Zeppelin, row::Int, dir::AbstractString="MAG", ext::AbstractString=".tif")

Returns the name of the spectrum/image file for the particle in the specified row.
"""
function spectrumfilename(zep::Zeppelin, row::Int, dir::AbstractString = "MAG", ext::AbstractString = ".tif", relocated=false)
    mag = hasproperty(zep.data, :MAG) ? convert(Int, trunc(zep.data[row, :MAG])) : 0
    # First check if a spectrum file exists
    tmp=repr(zep.data[row, :NUMBER])
    if relocated && isdir(joinpath(dirname(zep.headerfile),"RELOCATED"))
        for n in 6:-1:4
            fn = joinpath(dirname(zep.headerfile), "RELOCATED", repeat('0',max(0,n-length(tmp)))*tmp*ext)
            if isfile(fn)
                return fn
            end
        end
    end
    for n in 6:-1:4
        fn = joinpath(dirname(zep.headerfile), "$(dir)$(mag)", repeat('0',max(0,n-length(tmp)))*tmp*ext)
        if isfile(fn)
            return fn
        end
    end
    # Default case
    return joinpath(dirname(zep.headerfile), "$(dir)$(mag)", repeat('0',max(0,5-length(tmp)))*tmp*ext)
end

"""
    ParticleClassifier

A type that implements:

    classify(zep::Zeppelin, sr::ParticleClassifier)::CategoricalArray{String}
"""
abstract type ParticleClassifier end

"""
    spectrum(zep::Zeppelin, row::Int, withImgs = true)::Union{Spectrum, missing}

Returns the Spectrum (with images) associated with the particle at row.  If withImgs
is true, the associated image or images are read.
"""
function spectrum(zep::Zeppelin, row::Int, withImgs = true, relocated=true)::Union{Spectrum,Missing}
    file, at = spectrumfilename(zep, row, "MAG", ".tif", relocated), missing
    if isfile(file)
        try
            at = loadspectrum(ASPEXTIFF, file; withImgs = withImgs)
        catch err
            showerror(stderr, err)
            @info "$(file) does not appear to be a valid ASPEX spectrum TIFF."
        end
        try
            at[:BeamEnergy] = beamenergy(zep, get(at, :BeamEnergy, 20.0e3))
            at[:ProbeCurrent] = get(at, :ProbeCurrent, probecurrent(zep, 1.0))
            at[:Signature] = filter(kv->kv[2]>0.0, Dict(elm => zep.data[row, convert(Symbol,elm)] for elm in elms(zep)))
            at[:Name] = "$(name(zep))[$(zep.data[row, :NUMBER]), $(zep.data[row, :CLASS])]"
        catch
            @info "Error adding properties to ASPEX TIFF file."
        end
    end
    return at
end

function iszeppelin(filename::String)
    open(filename, enc"WINDOWS-1252") do ios
        seekstart(ios)
        if isequal(uppercase(String(read(ios,11))),"PARAMETERS=")
            readline(ios) # read the rest of the line
            return isequal(uppercase(readline(ios)),"HEADER_FMT=ZEPP_1")
        end
    end
    return false
end

function writeZep(zep::Zeppelin,  hdzfilename::String)
    remapcolumnnames = Dict(
        "NUMBER"=> "PART#\t1\tINT16",
        "FIELD" => "FIELD#\t1\tINT16",
        "MAGFIELD" => "MAGFIELD#\t1\tINT16",
        "XABS" => "X_ABS\tmm\tFLOAT",
        "YABS" => "Y_ABS\tmm\tFLOAT",
        "XDAC" => "X_DAC\t1\tINT16",
        "YDAC" => "Y_DAC\t1\tINT16",
        "XFERET" => "X_FERET\tµm\tFLOAT",
        "YFERET" => "Y_FERET\tµm\tFLOAT",
        "DAVG" => "DAVE\tµm\tFLOAT",
        "DMAX" => "DMAX\tµm\tFLOAT",
        "DMIN" => "DMIN\tµm\tFLOAT",
        "DPERP" => "DPERP\tµm\tFLOAT",
        "ASPECT" =>"ASPECT\t1\tFLOAT",
        "AREA" => "AREA\tµm²\tFLOAT",
        "PERIMETER" => "PERIMETER\tµm\tFLOAT",
        "ORIENTATION" => "ORIENTATION\tdeg\tFLOAT",
        "LIVETIME" => "LIVE_TIME\ts\tFLOAT",
        "FITQUAL" => "FIT_QUAL\t1\tFLOAT",
        "MAG" => "MAG\t1\tINT16",
        "VIDEO" => "VIDEO\t1\tINT16",
        "IMPORTANCE" => "IMPORTANCE\t1\tINT16",
        "COUNTS" => "COUNTS\t1\tFLOAT",
        "MAGINDEX" => "MAG_INDEX\t1\tINT16",
        "FIRSTELM" => "FIRST_ELEM\t1\tINT16",
        "SECONDELM" => "SECOND_ELEM\t1\tINT16",
        "THIRDELM" => "THIRD_ELEM\t1\tINT16",
        "FOURTHELM" => "FOURTH_ELEM\t1\tINT16",
        "COUNTS1" => "FIRST_CONC\tcounts\tFLOAT",
        "COUNTS2" => "SECOND_CONC\tcounts\tFLOAT",
        "COUNTS3" => "THIRD_CONC\tcounts\tFLOAT",
        "COUNTS4" => "FOURTH_CONC\tcounts\tFLOAT",
        "FIRSTPCT" => "FIRST_PCT\t%\tFLOAT",
        "SECONDPCT" => "SECOND_PCT\t%\tFLOAT",
        "THIRDPCT" => "THIRD_PCT\t%\tFLOAT",
        "FOURTHPCT" => "FOURTH_PCT\t%\tFLOAT",
        "TYPE4ET" => "TYPE(4ET#)\t1\tLONG",
        "VOIDAREA" => "VOID_AREA\tµm²\tFLOAT",
        "RMSVIDEO" => "RMS_VIDEO\t1\tINT16",
        "FITQUAL" => "FIT_QUAL\t1\tFLOAT",
        "VERIFIEDCLASS" => "VERIFIED_CLASS\t1\tINT16",
        "EDGEROUGHNESS" => "EDGE_ROUGHNESS\t1\tFLOAT",
        "COMPHASH" => "COMP_HASH\t1\tLONG",
        "CLASS" => "PSEM_CLASS\t1\tINT16",
        "TYPE_4ET_" => "Type[4ET]\t1\tLONG",
    )
    els = elms(zep)
    merge!(remapcolumnnames, Dict( uppercase(elm.symbol) => "$(uppercase(elm.symbol))\t%(k)\tFLOAT" for elm in els))
    merge!(remapcolumnnames, Dict( "U_$(uppercase(elm.symbol))_"  => "U[$(uppercase(elm.symbol))]\t%(k)\tFLOAT" for elm in els))
    headeritems = copy(zep.header)
    # add back the element tags
    for (i,elm) in enumerate(els)
        headeritems["ELEM$(i-1)"] = "$(elm.symbol) $(z(elm)) 1"
    end
    headeritems["ELEMENTS"] = "$(length(els))"
    if "CLASS" in names(zep.data)
        for (i, cn) in enumerate(zep.classnames)
            headeritems["CLASS$(i-1)"]=cn
        end
        headeritems["CLASSES"]= "$(length(zep.classnames))"
    end
    headeritems["TOTAL_PARTICLES"] = "$(size(zep.data,1))"
    # write out the header
    colnames = names(zep.data)
    open(hdzfilename, enc"WINDOWS-1252", "w") do ios
        println(ios, "PARAMETERS=$(length(headeritems)+size(zep.data,2)+2)")
        println(ios, "HEADER_FMT=ZEPP_1")
        foreach(hk->println(ios,"$(hk)=$(headeritems[hk])"), sort(collect(keys(headeritems))))
        println(ios, "PARTICLE_PARAMETERS=$(size(zep.data,2))")
        foreach(colnames) do col
            println(ios, get(remapcolumnnames, col, "$col\t1\t"*uppercase(repr(eltype(zep.data[:,col])))))
        end
    end
    pxzfilename = replace( hdzfilename, r".[h|H][d|D][z|Z]"=>".pxz")
    zd = DataFrame(zep.data)
    # Replace prettified items with the index (either into CLASS# in header or z(elm))
    for (ic, col) in enumerate(names(zd))
        if col in ("CLASS", "VERIFIEDCLASS", "FIRSTELM", "SECONDELM", "THIRDELM", "FOURTHELM" )
            coldata=zd[:,col]
            select!(zd ,Not(ic))
            if col in ( "CLASS", "VERIFIEDCLASS" )
                insertcols!(zd, ic, col => [ zc.index for zc in coldata ]) # CLASS
            elseif col in ("FIRSTELM", "SECONDELM", "THIRDELM", "FOURTHELM" )
                insertcols!(zd, ic, col => [ ze isa Element ? z(ze) : 0 for ze in coldata ]) # Element
            end
        end
    end
    CSV.write(pxzfilename, zd, delim="\t", missingstring="-", header=false)
end

function beamenergy(zep::Zeppelin, def=missing)
    number(v) = parse(Float64, match(r"([+-]?[0-9]+[.]?[0-9]*)",v)[1])
    val = def
    try
        v = get(zep.header, "ACCELERATING_VOLTAGE", missing)
        val = !ismissing(v) ? number(v)*1.0e3 : def
    catch
        # Ignore it...
    end
    return val
end

function probecurrent(zep::Zeppelin, def=missing)
    number(v) = parse(Float64, match(r"([+-]?[0-9]+[.]?[0-9]*)",v)[1])
    val = def
    try
        v = get(zep.header, "PROBE_CURRENT", missing)
        val = !ismissing(v) ? number(v) : def
    catch
        # Ignore it...
    end
    return val
end

function magdata(zep::Zeppelin, index::Int)
    h = split(zep.header["MAG_FMT"],isspace)
    v = split(zep.header["MAG$index"],isspace)
    return Dict( string(h[i])=>parse(Float64, v[i]) for i in 1:min(length(h),length(v)) )
end

"""
    randomsubset(zep::Zeppelin, maxRows::Int)

Creates an ordered random subselection of `maxRows` without replacement of the rows in the `zep` dataset.  If maxRows
is larger than the number of rows in `zep` then a UnitRange with all rows is returned.
"""
randomsubset(zep::Zeppelin, maxRows::Int) =
    return maxRows<size(zep.data,1) ? #
        sort(Random.shuffle(collect(eachparticle(zep)))[1:maxRows]) :
        eachparticle(zep)

"""
    residual(zep::Zeppelin, row::Int, withImgs=false)::Union{Spectrum,Missing}

Returns the residual spectrum for the particle at `row` or missing if it does not exist. Never returns images.
"""
function residual(zep::Zeppelin, row::Int, withImgs=false)::Union{Spectrum,Missing}
    res = missing
    try
        filename = spectrumfilename(zep, row, "Residual", ".msa")
        res = isfile(filename) ? loadspectrum(filename) : missing
    catch
        # Ignore errors, just return missing
    end
    return res
end

"""
    maxparticle(zep::Zeppelin, rows::Union{AbstractVector{Int},UnitRange{Int}})
    maxresidual(zep::Zeppelin, rows::Union{AbstractVector{Int},UnitRange{Int}})

Computes the maxparticle spectrum from the specified particles by `rows`.
"""
function maxparticle(zep::Zeppelin, rows::Union{AbstractVector{Int},UnitRange{Int}})
    mp, firstspec = missing, missing
    for spec in filter(s->!ismissing(s), map(r->spectrum(zep,r,false),rows))
        @assert spec isa Spectrum
        firstspec = ismissing(firstspec) ? spec : firstspec
        cx = NeXLSpectrum.counts(spec)
        mp = ismissing(mp) ? cx : max.(mp, cx)
    end
    # Now copy it out as a spectrum
    props = copy(firstspec.properties)
    props[:Name] = "MaxParticle"
    return Spectrum(firstspec.energy, mp, props)
end


"""
    maxresidual(zep::Zeppelin, rows::Union{AbstractVector{Int},UnitRange{Int}}=1:1000000)

Computes the max-residual spectrum from the specified particles by `rows`.  The residual must have been
previously computed in quantify(zep, ....)
"""
function maxresidual(zep::Zeppelin, rows::Union{AbstractVector{Int},UnitRange{Int}})
    mp, firstspec = missing, missing
    for spec in filter(s->!ismissing(s), map(r->residual(zep,r,false),rows))
        @assert spec isa Spectrum
        firstspec = ismissing(firstspec) ? spec : firstspec
        cx = NeXLSpectrum.counts(spec)
        mp = ismissing(mp) ? cx : max.(mp, cx)
    end
    # Now copy it out as a spectrum
    props = copy(firstspec.properties)
    props[:Name] = "MaxResidual"
    return Spectrum(firstspec.energy, mp, props)
end

"""
    rowsmax(zep::Zeppelin, col::Symbol; n::Int=10000000, classname::Union{Missing,AbstractString}=missing)

Returns the row indices associated with the `n` maximum values in the `col` column.

Examples:

    # Plot ten largest by :DAVG
    plot(zep, rowsmax(zep, :DAVG, 10))
    # data from 100 most Iron-rich sorted by DAVG
    rowsmax(DataFrame, zep, rowsmax(zep, :FE, 100), sortCol=:DAVG)
"""
function rowsmax(zep::Zeppelin, col::Symbol; n::Int=10000000, classname::Union{Missing,AbstractString}=missing)
    if ismissing(classname)
        d = collect(zip(1:250,zep.data[:,col]))
        s = sort(d, lt=(a1,a2)->isless(a1[2],a2[2]),rev=true)
        return getindex.(s[intersect(eachparticle(zep),1:n)],1)
    else
        cr = filter(r->zep.data[r,"CLASS"]==classname, eachparticle(zep))
        d = collect(zip(cr, zep.data[cr,col]))
        s = sort(d, lt=(a1,a2)->isless(a1[2],a2[2]),rev=true)
        return getindex.(s[intersect(eachindex(cr),1:n)],1)
    end
end
"""
    rowsmin(zep::Zeppelin, col::Symbol; n::Int=10000000, classname::Union{Missing,AbstractString}=missing)

Returns the row indices associated with the `n` minimum values in the `col` column.

Examples:

    # Plot ten smallest by :DAVG
    plot(zep, rowsmin(zep, :DAVG, 10))
    # data from 100 most Iron-rich sorted by DAVG
    rowsmin(DataFrame, zep, rowsmax(zep, :FE, 100), sortCol=:DAVG)
"""
function rowsmin(zep::Zeppelin, col::Symbol; n::Int=20, classname::Union{Missing,AbstractString}=missing)
    if ismissing(classname)
        d = collect(zip(eachparticle(zep),zep.data[:,col]))
        s = sort(d, lt=(a1,a2)->isless(a1[2],a2[2]),rev=false)
        return getindex.(s[intersect(eachparticle(zep),1:n)],1)
    else
        cr = filter(r->zep.data[r,"CLASS"]==classname, eachparticle(zep))
        d = collect(zip(cr, zep.data[cr,col]))
        s = sort(d, lt=(a1,a2)->isless(a1[2],a2[2]),rev=false)
        return getindex.(s[intersect(eachindex(cr),1:n)],1)
    end
end

"""
    rowsclass(zep::Zeppelin, classname::AbstractString; shuffle=false, n=1000000)

Returns the row indices associated with the specified `classname`.  If `shuffle` is true,
the row indices are shuffled so that a randomized n can be plucked off to plot or other.

Example:

    # Plot a randomized selection of 10 "Calcite" particles
    plot(zep, rowsclass(zep, "Calcite", shuffle=true, n=10), xmax=8.0e3)
"""
function rowsclass(zep::Zeppelin, classname::AbstractString; shuffle=false, n=1000000)
    res = filter(i->zep.data[i,"CLASS"]==classname, eachparticle(zep))
    return (shuffle ? Random.shuffle(res) : res)[intersect(eachindex(res),1:n)]
end

rowsclass(zep::Zeppelin, classnames::AbstractVector{<:AbstractString}; shuffle=false, n=1000000) =
    mapreduce(cn->rowsclass(zep, cn, shuffle=shuffle, n=n), union, classnames)



"""
    Base.filter(filt::Function, zep::Zeppelin)

Use a function of the form `filt(row)::Bool` to filter a Zeppelin dataset returning a new Zeppelin dataset
with only the rows for which the function evaluated true.


Example:

    gsr = filter(row->startswith( String(row["CLASS"]), "GSR."), zep)
"""
function Base.filter(filt::Function, zep::Zeppelin)::Zeppelin 
    zd = filter(filt, zep.data; view=false)
    return Zeppelin(zep.headerfile, copy(zep.header), zd, copy(zep.classnames))
end

const MORPH_COLS = ( "NUMBER", "XABS", "YABS", "DAVG", "DMIN", "DMAX", "DPERP", "PERIMETER", "ORIENTATION", "AREA" )
const CLASS_COLS = ( "CLASS", "VERIFIEDCLASS", "IMPORTANCE" )
const COMP_COLS = ( "FIRSTELM", "FIRSTPCT", "SECONDELM", "SECONDPCT", "THIRDELM", "THIRDPCT", "FOURTHELM", "FOURTHPCT", "COUNTS" )
const ALL_ELMS = map(elm->uppercase(elm.symbol), elements[1:95])

const ALL_COMPOSITIONAL_COLUMNS = ( "FIRST", "FIRSTELM", "SECONDELM", "THIRDELM", "FOURTHELM", "FIRSTELM", "SECONDELM",
    "THIRDELM", "FOURTHELM", "COUNTS1", "COUNTS2", "COUNTS3", "COUNTS4", "FIRSTPCT", "SECONDPCT", "THIRDPCT",
    "FOURTHPCT", "TYPE4ET", "COUNTS", "FITQUAL", "COMPHASH" )
const ALL_CLASS_COLS = ( "CLASS", "VERIFIEDCLASS", "IMPORTANCE" )


"""
    translate(zep::Zeppelin, am::AbstractAffineMap)::Zeppelin

The result is a copy of `zep` with the columns `:XABS` and `:YABS` replaced with the DataFtranslated coordinates.
"""
function translate(zep::Zeppelin, am::AbstractAffineMap)::Zeppelin
    res=copy(zep)
    foreach(eachrow(res.data)) do r
        ( r.XABS, r.YABS ) = am(SVector(r.XABS, r.YABS))
    end
    return res
end

"""
    align(zep1::Zeppelin, zep2::Zeppelin; tol=0.001, finealign=true)

Aligns `zep2` to overlay with `zep1`.   `zep1` and `zep2` are assumed to be particle data sets collected
from the same sample but that may translated and rotated from one-another.   Not all the particles need to
correspond between `zep1` and `zep2` but many must.
The function returns a copy of `zep2` transformed to overlay `zep1`.
"""
function align(zep1::Zeppelin, zep2::Zeppelin; tol=0.001, finealign=true)
    ps1 = map(xy-> SA[ xy... ], zip(zep1[:, :XABS], zep1[:,:YABS]))
    ps2 = map(xy-> SA[ xy... ], zip(zep2[:, :XABS], zep2[:,:YABS]))
    ct1, ct2 = align(ps1, ps2, tol=tol, finealign=finealign)
    return translate(zep2, inv(ct1)∘ct2)
end
"""
    correspondences(zep1::Zeppelin, zep2::Zeppelin; tol=0.01, invert=false)

Identify particle correspondences between `zep1` and `zep2`.  Returns two new
`Zeppelin` data sets `zc1` and `zc2` in which corresponding particles are matched 
by row index. So `zc1[1,:]` is likely to refer to the same particle as `zc2[1, :]`.
Unless invert is chosen, in which case only the particles with no corresponding 
partners are returned.
"""
function correspondences(zep1::Zeppelin, zep2::Zeppelin; tol=0.01, invert=false)
    ps1 = map(xy-> SA[ xy... ], zip(zep1[:, :XABS], zep1[:,:YABS]))
    ps2 = map(xy-> SA[ xy... ], zip(zep2[:, :XABS], zep2[:,:YABS]))
    c1, c2 = correspondences(ps1, ps2, tol=tol, invert=invert)
    return ( 
        Zeppelin(zep1.headerfile, copy(zep1.header), zep1.data[c1, :], copy(zep1.classnames)),
        Zeppelin(zep2.headerfile, copy(zep2.header), zep2.data[c2, :], copy(zep2.classnames))
    )
end


"""
    identify(zeps::AbstractArray{Zeppelin}; tol=0.001, ctol=0.01, columns=())::DataFrame

Takes multiple Zeppelin data sets that represent the same particles and returns a `DataFrame`
that identifies the particles from one data set to the next.  The algorithm aligns the data sets
by determining the rotation and offset to bring them into registration.  Then it identifies all 
the unique particles by position and tracks them from one data set to the next.  The number of
rows in the `DataFrame` is the number of unique particle positions identified (the `:XABS` and
`:YABS` colums).  The `:COUNT` column is the number of times a particle was found at this 
position.  The `:APPEARS` column is the first particle data set in which it was found by index.
The first columns are the index at which the particle is found in the `Zeppelin` data set.

The columns argument allows you to extract `mean` and `std` for a property of a physical particle 
within the data sets in which it was measured.  Like: `columns=(:DAVG, :AREA)`
"""
function identify(zeps::AbstractArray{Zeppelin}; tol=0.001, ctol=0.01, columns=())
    pss=map(zeps) do zep
        map(xy-> SA[ xy... ], zip(zep[:, :XABS], zep[:,:YABS]))
    end
    res = identify(pss; tol=tol, ctol=ctol)
    rename!(res, ("PS$i"=>NeXLSpectrum.name(zeps[i]) for i in eachindex(zeps))...)
    l = levels(res[:, :APPEARS])
    insertcols!(res, length(zeps)+2, :ZEPPELIN => map(i->name(zeps[l[i.ref]]), res[:,:APPEARS]))
    for col in columns
        cm=map(eachrow(res)) do r
            mean(skipmissing( [ismissing(r[i]) ? missing : z[r[i], col] for (i, z) in enumerate(zeps)]))
        end
        cs=map(eachrow(res)) do r
            std(skipmissing( [ismissing(r[i]) ? missing : z[r[i], col] for (i, z) in enumerate(zeps)]))
        end
        insertcols!(res, Symbol(col,"_mean")=>cm)
        insertcols!(res, Symbol(col,"_std")=>cs)
    end
    res
end

Base.eachrow(zep::Zeppelin) = eachrow(zep.data)
DataAPI.nrow(zep::Zeppelin) = nrow(zep.data)
DataAPI.ncol(zep::Zeppelin) = ncol(zep.data)
function Base.sort!(zep::Zeppelin, cols=All(); nargs...)
    sort!(zep.data, cols; nargs...)
    return zep
end
function Base.sort(zep::Zeppelin, cols=All(); nargs...)
    sd = sort(zep.data, cols; nargs...)
    return Zeppelin(zep.headerfile, copy(zep.header), sd, copy(zep.classnames))
end

"""
   reclass!(z::Zeppelin, r::Int, new_name::String)
   
Assigns the particle in row `r` to a class named `new_name`.
If the class name isn't already available, it is added to the
list of available classes.
"""
function reclass!(z::Zeppelin, r::Int, new_name)
    nn = String(new_name)
    i = findfirst(n->n==nn, z.classnames)
    if isnothing(i)
        push!(z.classnames,nn)
        i = length(z.classnames)
    end
    nc = NeXLParticle.ZepClass(z,i-1)
    z.data[r,:CLASS] = nc
    nc
end


"""
   reclassall!(z1::Zeppelin, z2::Zeppelin, old_name)

Re-assign the class for all rows with CLASS=old_name to
the CLASS of the same particle in `z2`.
"""
function reclassall!(z1::Zeppelin, z2::Zeppelin, old_name)
    foreach(eachparticle(z1)) do i
        if z1[i,:CLASS]==old_name
            reclass!(z1,i,z2[i,:CLASS])
        end
    end
end
