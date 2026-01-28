# Open a PNG graphics device
png("test_plot.png", width = 800, height = 600)

# Make a simple plot
x <- 1:10
y <- x^2

plot(x, y,
     type = "b",
     col = "blue",
     pch = 19,
     main = "Conda R Environment Test",
     xlab = "X values",
     ylab = "X squared")

# Close the graphics device (this actually writes the file)
dev.off()

cat("Plot successfully created: test_plot.png\n")
