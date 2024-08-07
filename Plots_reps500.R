library(png)

# Initialize an empty list to store filenames
filenames <- list()

# Define the parameters
eps <- 0:5
timesteps <- c("1", "01", "001", "0001")
reps <- 500

# Generate the filenames
for (ep in eps) {
  for (dt in timesteps) {
      fname <- paste("eps", ep, "dt", dt, "reps", reps, ".png", sep = "")
      filenames <- append(filenames, fname)
    }
}

# Print the filenames
print(filenames)

images <- lapply(filenames, function(file) {
  if (file.exists(file)) {
    readPNG(file)
  } else {
    # If the file doesn't exist, return NULL
    print(file)
    NULL
  }
})

# Remove NULL entries (files that don't exist)
images <- Filter(Negate(is.null), images)

# Ensure we have exactly 24 images for a 4x6 grid
#if (length(images) != 24) {
  #stop("The number of images is not 24. Adjust the filenames or parameters.")
#}


# Define the output file name and dimensions
output_file <- "compare_image.png"
width <- 2400      # Width of the output image in pixels
height <- 2000     # Height of the output image in pixels
res <- 600         # Resolution in DPI (dots per inch)
# Set up the png device with specified dimensions and resolution
png(filename = output_file, width = width, height = height, res = res)

# Set up the plotting area for a 4x6 grid
par(mfrow = c(6,4), mar = c(0, 0, 0, 0))

# Plot each image in the grid
for (img in images) {
  plot.new()
  rasterImage(img, 0, 0, 1, 1)
}
dev.off()

