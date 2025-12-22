# Number of iterations
total_iterations <- 100

# Create a progress bar
pb <- utils::txtProgressBar(min = 0, max = 2, style = 3)

# Simulate a process with a for loop
for (i in 1:total_iterations) {
  Sys.sleep(0.1)  # Simulate a delay for each iteration
  utils::setTxtProgressBar(pb, i)  # Update the progress bar
}

# Close the progress bar when done
close(pb)
