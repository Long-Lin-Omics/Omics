# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_file <- args[2]
sample_name <- args[3]

# Read fragment lengths from input file
data <- read.table(input_file, header=FALSE)

# Calculate histogram data
#hist_data <- hist(data$V1, breaks=10, plot=FALSE)

# Plot as a line plot and save to specified output file

png(output_file, width=1200, height=600)
par(mfrow=c(1,2))

# Plot the full range histogram
hist(data$V1, breaks=10, col="blue", 
     xlab="Fragment Length (bp)", ylab="Frequency",
     main=paste0("Fragment Length Distribution of ",sample_name))

# Plot the zoomed-in version with x-axis limited to 2000 bp
hist(data$V1[data$V1 <= 2000], col="blue", 
     xlab="Fragment Length (bp)", ylab="Frequency",
     main=paste0("Fragment Length Distribution [0, 2k] of ",sample_name))
     
dev.off()