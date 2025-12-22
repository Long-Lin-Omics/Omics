library(hdf5r)

# Function to extract the read by ID from an input FAST5 and write it to a new FAST5
extract_read_to_file <- function(input_fast5, output_fast5, read_id) {
  # Open the input FAST5 file
  file_in <- H5File$new(input_fast5, mode = "r")
  
  # Check if the read_id exists in the file
  if (read_id %in% names(file_in)) {
    # Open the output FAST5 file for writing
    file_out <- H5File$new(output_fast5, mode = "w")
    
    # Copy the read group to the new file
    tree = file_in$ls(recursive=T)
    tree_read = tree[grepl(read_id,tree$name),]
    for (i in seq(1,dim(tree_read)[1])){
      name = tree_read[i, 'name']
      cat(name,'\n')
      # H5I_DATASET is not a group
      if (tree_read[i, 'obj_type'] == 'H5I_DATASET'){
        file_out[[name]] = file_in[[name]]$read()
      }
      else{
        file_out$create_group(name)  
      }
      # every node can have attrs no matter groups or datasets
      if (tree_read[i, 'num_attrs'] !=0 ){
       old_attr_names = h5attr_names(file_in[[name]])
       for (attr in old_attr_names){
        h5attr(file_out[[name]],attr) = h5attr(file_in[[name]],attr)
       }
      }
      
    }
    
    cat(sprintf("Read %s successfully extracted to %s\n", read_id, output_fast5))
    
    # Close the output file
    file_out$close()
  } else {
    cat(sprintf("Read with ID %s not found in %s\n", read_id, input_fast5))
  }
  
  # Close the input file
  file_in$close()
}

# Example usage:
input_fast5 <- "/ddn/gs1/project/nextgen/post/hug4/Nanopore_Data/projects/E14_PCR_Cycle_TEST/basecalls4/workspace/fast5//FAO78096_8ed9fe38_0.fast5"
output_fast5 <- "path_to_output.fast5"
read_id <- "read_638e2a2e-b710-4c58-8bf8-33d7772f2341"  # Replace with the read ID you are interested in

# Check if the script is being run from the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 3) {
  # Command line parameters are provided
  input_fast5 <- args[1]
  output_fast5 <- args[2]
  read_id <- args[3]
  if (!startsWith(read_id,'read_')){
    read_id = paste0('read_',read_id)
  }
}

extract_read_to_file(input_fast5, output_fast5, read_id)
