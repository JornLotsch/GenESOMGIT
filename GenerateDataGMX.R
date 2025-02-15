GeneratedDataGMX <- function(Data, DensityRadius, Cls, GenPerData = 10, PlotIt = 0, Plot3D = FALSE, AlphaGenerated = 0.4) {
  # Generate New data with the same structure of Data given a valid U- and P-matrix
  #
  # INPUT
  # Data(1:n,1:d)             data matrix, cases in rows, features in columns
  # DensityRadius             Radius > 0 used for generating a suitable P-Matrix
  #
  # OPTIONAL
  # GenPerData                      number of newly generated data points per given data point
  # dafault  GenPerData =10;
  # PlotIt                          = 1 Plot the sigmoid, = 2: Plot data density, defaut: PlotIt = 0;
  # OUTPUT
  # GeneratedData(1:GenPerData*n,d) the generated data
  # Umatrix(1:Lines,1:Columns)     a matrix of U-heigts such that Data is located at BestMatches positions.
  # Wts(1:Lines*Columns,1:d)       weight vectors of neurons in retina odered such that savewts(... Wts,...) produces
  #                                a valid *.wts file and is consistent with the other outputs

  # AUTHOR: ALU 2023, R version: JL 2023

  # Parameter for data generation
  MaxD <- 2.0 # Data will be generated up to r*MaxD
  LimitAB <- 0.72 # AB limit in the sigmoid
  LimitBC <- 1.22 # BC limit in the sigmoid
  PercentA <- 0.80 # fraction of data in set A  (core)
  PercentB <- 0.15 # fraction of data in set  B(decreasing)
  PercentC <- 1 - (PercentA + PercentB)

  r <- 0.4 * DensityRadius # shrink radius to below half of the intercluster distance

  n <- dim(Data)[1] # n = nr of data Points
  d <- dim(Data)[2] # d = dimensionality of data

  # Class labeling. Check if provided, else all cases in class 1
  if (hasArg("Cls") == TRUE) {
    if (length(Cls) != n) {
      stop("Unequal number of cases and class labels.")
    }
  } else {
    Cls <- rep(1, n)
  }


  N <- n * GenPerData # number of generated data set instances
  Jitter <- matrix(rnorm(N * d, mean = 0, sd = 1.5), N, d) * r #   generate hypercube
  L <- sqrt(rowSums(Jitter^2)) # length of vectors
  Ind <- which(L > 0)
  Jitter[Ind, ] <- Jitter[Ind, ] / (L[Ind]) # normalize to unified sphere => all have length = 1

  AnzA <- round(N * PercentA)
  AnzB <- round(N * PercentB)
  AnzC <- round(N * PercentC)
  IndA <- 1:AnzA
  IndB <- (AnzA + 1):(AnzA + AnzB)
  IndC <- (AnzA + AnzB + 1):(AnzA + AnzB + AnzC) # IndC = min(N,IndC);
  AnzA <- length(IndA)
  AnzB <- length(IndB)
  AnzC <- length(IndC)
  AnzABC <- AnzA + AnzB + AnzC

  # Generate lengths of vectors for sets A, B and C
  SigmoidLength <- c(runif(AnzA, 0, LimitAB), runif(AnzB, LimitAB, LimitBC), runif(AnzC, LimitBC, MaxD))

  SigmoidLength <- SigmoidLength * r # rescale from 1 to r

  # Extend vectors
  SigmoJid <- matrix(rep(SigmoidLength, d), ncol = d)

  Jidder <- Jitter * SigmoJid
  GenInd <- rep(1:n, GenPerData) # GenPerdata mal [1:n]
  GeneratedData <- Data[GenInd, ] + Jidder

  GenCls <- rep(Cls, times = GenPerData) # Classes of generated data

  # Plot data
  pSigmoid <- NULL
  pPDEsigmoid <- NULL
  pGeneratedData <- NULL
  if (PlotIt > 0) { # plot sigmoid and data distribution

    source("/home/joern/Dokumente/Aktuell/PDEplotGG.R")

    ProcA <- round(AnzA / AnzABC * 100)
    ProcB <- round(AnzB / AnzABC * 100)
    ProcC <- round(AnzC / AnzABC * 100, 1)

    # Parameter of the sigmoid
    Dx <- seq(0, 2, 0.01)
    Tx <- 0.1
    Sigmoid <- 1 - 1. / (1 + exp(-(Dx - 1) / Tx))

    require(ggplot2)
    y <- 0.15
    pSigmoid <-
      ggplot(data = cbind.data.frame(D = Dx, Sigmoid = Sigmoid), aes(x = D, y = Sigmoid)) +
      geom_line(color = "dodgerblue") +
      geom_vline(xintercept = c(LimitAB, LimitBC), color = "darkgreen") +
      geom_vline(xintercept = MaxD, color = "salmon") +
      theme_light() +
      annotate(
        geom = "text", x = c(0.32, 0.94, 1.52), y = y, label = c(paste0(LETTERS[1:3], " ", c(ProcA, ProcB, ProcC), "#")),
        color = "red"
      ) +
      annotate(
        geom = "text", x = 0, y = 0.9, label = paste0("d = 1 at DensityRadius r = ", DensityRadius),
        color = "black", hjust = 0
      ) +
      labs(
        title = "prob(genData is neighbor of data point)",
        x = "Distance / sphere radius"
      )

    pPDEsigmoid <-
      PDEplotGG(SigmoidLength) +
      geom_vline(xintercept = c(LimitAB * r, LimitBC * r), color = "darkgreen") +
      geom_vline(xintercept = MaxD * r, color = "salmon") +
      theme_light() +
      labs(title = "prob(Generated Data) vs Dist to Data Point", x = "Distances", y = "pdf") +
      guides(color = "none")
  }

  if (PlotIt > 1) { # plot data density
    if (d >= 1 & Plot3D == FALSE) {
      XY <- cbind.data.frame(X = GeneratedData[, 1], Y = GeneratedData[, 2], GenCls = GenCls)
      XYorig <- cbind.data.frame(X = Data[, 1], Y = Data[, 2], Cls = Cls)
      pGeneratedData <-
        ggplot() +
        geom_point(data = XY, aes(x = X, y = Y, color = factor(GenCls)), alpha = AlphaGenerated) +
        geom_point(data = XYorig, aes(x = X, y = Y, color = factor(Cls))) +
        theme_light() +
        labs(title = "Data density", color = "Class")
    } else {
      if (d > 2) {
        require(scatterplot3d)
        require(gridGraphics)
        require(grid)
        require(gridExtra)

        grab_grob <- function() {
          grid.echo()
          grid.grab()
        }

        XYZ <- cbind.data.frame(X = GeneratedData[, 1], Y = GeneratedData[, 2], Z = GeneratedData[, 3], GenCls = factor(GenCls))
        scatterplot3d(x = XYZ$X, y = XYZ$Y, z = XYZ$Z, color = as.integer(as.factor(XYZ$GenCls)), main = "Data density")
        pGeneratedDensity <- grab_grob()
      }
    }
  }
  return(list(
    OriginalData = Data,
    OriginalClasses = Cls,
    GeneratedData = GeneratedData,
    GeneratedClasses = GenCls,
    pSigmoid = pSigmoid,
    pPDEsigmoid = pPDEsigmoid,
    pGeneratedData = pGeneratedData
  ))
}
