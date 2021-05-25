using VideoIO
using ProgressMeter
using Plots

# Load image and create the video object
io = VideoIO.open("../../../../Downloads/video.mp4")
f = VideoIO.openvideo(io)

# Create the necessary objects to read frames and show progress
img = read(f)
current = Set(img)
cache_size = Int64[]
cache_throughput = Int64[]
frames = VideoIO.get_number_frames("../../../../Downloads/video.mp4")
p = Progress(frames, 1, "Analysing video", 50)
next!(p)

# Read the file until the end of the video
while !eof(f)
    read!(f, img)
    new = Set(img)
    append!(cache_size, length(new))
    append!(cache_throughput, length(setdiff(new, current)))
    current = copy(new)
    next!(p)
end
close(f)

# Plot the required cache size and print its mean
Plots.plot(cache_size * 3)
Plots.xlabel!("Frame")
Plots.ylabel!("Cache size in bytes", fontfamily="Computer Modern", legend=false)
Plots.savefig("../../final_report/figures/cache_size.pdf")
mean(cache_size * 3)

# Plot the the required cache throughput and print its mean
Plots.plot(cache_throughput .* 3 ./ 25)
Plots.xlabel!("Frame")
Plots.ylabel!("Cache throughput in bytes/s", fontfamily="Computer Modern", legend=false)
Plots.savefig("../../final_report/figures/cache_throughput.pdf")
mean(cache_throughput) .* 3 ./ 25
