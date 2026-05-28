
using Dates
using CSV
using DataFrames

loaded_data = DataFrame()
loaded_time = [DateTime(0), DateTime(0)]

function get_frequency(t_start::DateTime, period::Real)
    t_end = t_start + Second(period)
    
    if isempty(loaded_data) || t_start < loaded_time[1] || t_end > loaded_time[2]
        load_period(t_start, t_end)
        @show size(loaded_data)

    end

    # Filter the loaded data for the requested period
    df = loaded_data[(loaded_data.t .>= t_start) .& (loaded_data.t .<= t_end), :]

    return df
end


function load_period(t_start::DateTime, t_end::DateTime)
    if t_end < t_start
        error("End time must be after start time")
    end
    if t_start < Date(2014, 10, 1) || t_end > Date(2024, 11, 30)
        error("GridFrequency: Data for date not available")
    end
    
    start_date = Dates.format(t_start,"yyyy_mm")
    filename = joinpath(@__DIR__, "RTE_Frequence_$(start_date).txt")

    start_month = CSV.read(filename, DataFrame, header=[:t, :f], skipto=2, delim=';', dateformat="dd/mm/yyyy HH:MM:SS", ignorerepeated=true)
    global loaded_data
    global loaded_time

    if Dates.month(t_start) == Dates.month(t_end)
        loaded_data = start_month
        @show size(loaded_data)
    else
        end_date = Dates.format(t_end,"yyyy_mm")
        filename = joinpath(@__DIR__, "RTE_Frequence_$(end_date).txt")
        end_month = CSV.read(filename, DataFrame, header=[:t, :f], skipto=2, delim=';', dateformat="dd/mm/yyyy HH:MM:SS", ignorerepeated=true)
        loaded_data = vcat(start_month, end_month)
    end

    loaded_time = [minimum(loaded_data.t), maximum(loaded_data.t)]
end
