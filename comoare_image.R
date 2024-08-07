library(png)
library(jpeg)  # For reading JPG files

# Initialize an empty list to store filenames
filenames <- list()

# Define the parameters
eps <- 0:5
timesteps <- c("1", "01", "001", "0001")
reps <- 1

# Generate the filenames and handle .JPG only if .png does not exist
for (ep in eps) {
  for (dt in timesteps) {
    # Generate .png and .JPG filenames
    png_fname <- paste("eps", ep, "dt", dt, "reps", reps, ".png", sep = "")
    jpg_fname <- paste("eps", ep, "dt", dt, "reps", reps, ".JPG", sep = "")
    
    # Check for .png file first
    if (file.exists(png_fname)) {
      filenames <- append(filenames, png_fname)
    } else if (file.exists(jpg_fname)) {
      filenames <- append(filenames, jpg_fname)
    } else {
      # Print missing filenames for debugging
      print(paste("File not found:", png_fname, jpg_fname))
    }
  }
}

# Print the filenames for verification
print(filenames)

# Read all images into a list
images <- lapply(filenames, function(file) {
  if (file.exists(file)) {
    if (grepl("\\.png$", file)) {
      readPNG(file)
    } else if (grepl("\\.JPG$", file)) {
      print(file)
      readJPEG(file)
    }
  } else {
    # If the file doesn't exist, return NULL
    print(file)
    NULL
  }
})

# Remove NULL entries (files that don't exist or could not be read)
images <- Filter(Negate(is.null), images)

# Check the number of images loaded
if (length(images) != 24) {
  stop("The number of images is not 24. Adjust the filenames or parameters.")
}

# Define the output file name and dimensions
output_file <- "final_image.png"
width <- 2400  # Width of the output image in pixels
height <- 2400 # Height of the output image in pixels
res <- 600     # Resolution in DPI

# Set up the png device with specified dimensions and resolution
png(filename = output_file, width = width, height = height, res = res)

# Set up the plotting area for a 4x6 grid
par(mfrow = c(6, 4), mar = c(0, 0, 0, 0))

# Plot each image in the grid
for (img in images) {
  plot.new()
  rasterImage(img, 0, 0, 1, 1)
}

# Turn off the device
dev.off()
