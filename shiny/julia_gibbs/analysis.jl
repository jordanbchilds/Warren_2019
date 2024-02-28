using Pkg;

using Statistics, Random, NamedArrays, DataStructures, CSV, DataFrames, GLM, StatsBase, Distributions;

include("./functions.jl")

fullraw_data = CSV.read("../dat.txt", DataFrame)
raw_data = fullraw_data[:,["value", "id", "channel", "patient_type", "patient_id", "cell_id"]];

mitochan = "VDAC1"
channels = ["NDUFB8", "NDUFB13", "SDHA", "UqCRC2", "COX4+4L2", "MTCO1", "OSCP"]
nChains = length(channels) ;

raw_data[ raw_data[:,"channel"].=="GRIM19", "channel"] .= "NDUFB13";

data_ind = [ xx in [mitochan; channels] for xx in raw_data[:,"channel"]]
data = raw_data[data_ind, :]
data[!,"value"] = log.(data[!,"value"] ) ;

sbjIDs = unique(data[:,"patient_id"])
ctrlIDs = unique( data[ data[:,"patient_type"].=="control", "patient_id"] )
ptsIDs = unique( data[ data[:,"patient_type"].=="patient", "patient_id"] ) |> sort
ptsIDs_ind = [ !(xx in ["P01", "P02"]) for xx in ptsIDs]
ptsIDs = ptsIDs[ ptsIDs_ind ];

rename!(data, ["Value", "fibreID", "Channel", "sbjType", "sampleID", "cellID"]) ;

if !ispath("Output")
    mkdir("Output")
end

nChains = 5

Threads.@threads for chan in channels
    Threads.@threads for pat in ptsIDs
        root = chan*"__"*pat
        dd = getData_mats(data; mitochan="VDAC1", chan=chan, pts=[pat], ctrlIDs=ctrlIDs)
        Threads.@threads for i in 1:nChains
           output::Dict{String, NamedArray} = gibbs_sampler(dd, warmup=100, iter=100)
           mySaver(output, fileRoot=string("Output/"*root *"_chain"*lpad(i, 2, "0")*"_") )
        end
    end
end