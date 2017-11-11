# EEMindices: Automates the calculation of indices derived from excitation-emission matrix (EEM) spectra of dissolved organic matter (DOM).
#
# Copyright © 2017 Suzanne Hodgkins and Florida State University.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# For details on the use of this program, see the included README file.
#
# Contact: Dr. Suzanne Hodgkins, suzanne.b.hodgkins@gmail.com

# Startup notice ----
cat("EEMindices: Copyright © 2017 Suzanne Hodgkins and Florida State University.
This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you
are welcome to redistribute it under the terms of the GNU General Public
License (version 3 or later). This license is included with this program and
is also available at <http://www.gnu.org/licenses/>.")

# Function ----
library(colorRamps)
EEMindices <- function(filepattern, outputfile='EEMS_results.csv', show.plots=TRUE, fixed=TRUE) {

  if(fixed) {
    eem.files <- list.files()
    eem.files <- eem.files[grep(filepattern, eem.files, fixed=TRUE)]
  } else {
    eem.files <- list.files(pattern=filepattern)
  }

  if(length(grep('\\.txt', eem.files))==0) {
    stop('Files should be saved in .txt format.')
  }

  sample.names <- gsub(filepattern, "", eem.files)
  sample.names <- gsub("\\.txt", "", sample.names)

  # Preallocate results arrays:
  eem.matrices <- setNames(as.list(rep(0, length(eem.files))), eem.files)
  BIX <- rep(0, length(eem.files))
  FI <- rep(0, length(eem.files))
  HIX <- rep(0, length(eem.files))
  C.intensity <- rep(0, length(eem.files))
  C.em.max <- rep(0, length(eem.files))
  C.ex.max <- rep(0, length(eem.files))
  ini.C.em.max <- rep(0, length(eem.files))
  min.C.ex <- rep(0, length(eem.files))
  testcount <- rep(0, length(eem.files))
  notes <- rep(0, length(eem.files))

  trapz <- function(x, y) { 2 * (y[1]/2 + sum(y[2:(length(y)-1)]) + y[length(y)]/2) }
  # trapezoid-rule integration function for calculating HIX in EEM files

  for(i in seq_along(eem.matrices)) {
    eem.matrices[[i]] <- data.matrix(read.table(file=names(eem.matrices)[i], header=TRUE, skip=12))
    dimnames(eem.matrices[[i]])[[2]] <- gsub("X", "", dimnames(eem.matrices[[i]])[[2]])  # get rid of 'X' in column names
    dimnames(eem.matrices[[i]])[[2]] <- gsub("\\.0", "", dimnames(eem.matrices[[i]])[[2]])  # get rid of '.0' in column names
    dimnames(eem.matrices[[i]])[[1]] <- gsub("\\.0", "", dimnames(eem.matrices[[i]])[[1]])  # get rid of '.0' in row names

    em380ex310 <- eem.matrices[[i]]['380', '310']
    em430ex310 <- eem.matrices[[i]]['430', '310']
    BIX[i] <- em380ex310/em430ex310

    em470ex370 <- eem.matrices[[i]]['470', '370']
    em520ex370 <- eem.matrices[[i]]['520', '370']
    FI[i] <- em470ex370/em520ex370

    H <- trapz(seq(434, 480, by=2), eem.matrices[[i]][as.character(seq(434, 480, by=2)), '255'])
    L <- trapz(seq(300, 344, by=2), eem.matrices[[i]][as.character(seq(300, 344, by=2)), '255'])
    HIX[i] <- H/L

    # Optimize window for finding the location of the C peak maximum
    # (1) Define an initial emission wavelength, with high C fluorescence intensity, for finding the trough between the A and C peaks:
    ini.C.em.max[i] <- dimnames(eem.matrices[[i]])[[1]][which(eem.matrices[[i]][,'325']==max(eem.matrices[[i]][as.character(seq(376, 524, by=2)),'325']))]
    # (2) Find the excitaton wavelength of the trough between the A and C peaks:
    range.trough <- eem.matrices[[i]][ini.C.em.max[i],9:18]
    is.local.min <- (range.trough < eem.matrices[[i]][ini.C.em.max[i],10:19]) & (range.trough < eem.matrices[[i]][ini.C.em.max[i],8:17])
    if(any(is.local.min)) {
      min.C.ex[i] <- names(which(range.trough[is.local.min]==min(range.trough[is.local.min])))
    } else {
      #min.C.ex[i] <- names(which(range.trough==min(range.trough))) # old definition
      min.C.ex[i] <- '300' # new definition (may adjust this value)
      warning(sprintf("Failed to find the A-C boundary for sample:  %s", sample.names[i]))
    }  # minimum (must be local minimum, if possible) between excitation wavelengths 250-325nm at emission=ini.C.em.max
    rm(range.trough, is.local.min)

    min.C.index <- (as.numeric(min.C.ex[i]) - 235)/5  # array index corresponding to excitation wavelength of A-C boundary trough
    C.em.range.index <- ((as.numeric(ini.C.em.max[i])-26-288)/2):((as.numeric(ini.C.em.max[i])+26-288)/2) # array indices corresponding to approximate range of C emission wavelengths

    # Find the C peak intensity maximum, and the corresponding emission and excitation wavelengths, and check the results:
    testable.range <- eem.matrices[[i]][C.em.range.index, min.C.index:33] # check for C.intensity within range determined above
    testable.values <- as.vector(testable.range)
    testcount[i] <- 1          # keep track of number of tests

    while(testcount[i]<=100) {   # Do no more than 100 tests.

      C.intensity[i] <- max(testable.values)
      wavelength.C.em.max <- dimnames(testable.range)[[1]][which(testable.range==C.intensity[i], arr.ind=TRUE)[1,1]]
      wavelength.C.ex.max <- dimnames(testable.range)[[2]][which(testable.range==C.intensity[i], arr.ind=TRUE)[1,2]]

      if(nrow(which(testable.range==C.intensity[i], arr.ind=TRUE)) > 1) {
        warning(paste0(sample.names[i], ": Multiple cells in the EEM equal the C intensity."))
      }
      index.C.em.max <- which(dimnames(eem.matrices[[i]])[[1]]==wavelength.C.em.max)
      index.C.ex.max <- which(dimnames(eem.matrices[[i]])[[2]]==wavelength.C.ex.max)

      # if(nrow(which(eem.matrices[[i]]==C.intensity[i], arr.ind=TRUE)) > 1) {
      #   warning(paste0(sample.names[i], ": Multiple cells in the EEM equal the C intensity."))
      # }
      # index.C.em.max <- which(eem.matrices[[i]]==C.intensity[i], arr.ind=TRUE)[1,1]
      # index.C.ex.max <- which(eem.matrices[[i]]==C.intensity[i], arr.ind=TRUE)[1,2]

      C.em.max[i] <- dimnames(eem.matrices[[i]])[[1]][index.C.em.max]
      C.ex.max[i] <- dimnames(eem.matrices[[i]])[[2]][index.C.ex.max]

      # make sure the calculated C.intensity is actually higher than all of the neighboring intensities
      max.neighbor.intensity <- max(eem.matrices[[i]][index.C.em.max, index.C.ex.max],
                                    eem.matrices[[i]][index.C.em.max, index.C.ex.max+1],
                                    eem.matrices[[i]][index.C.em.max, index.C.ex.max-1],
                                    eem.matrices[[i]][index.C.em.max+1, index.C.ex.max],
                                    eem.matrices[[i]][index.C.em.max+1, index.C.ex.max+1],
                                    eem.matrices[[i]][index.C.em.max+1, index.C.ex.max-1],
                                    eem.matrices[[i]][index.C.em.max-1, index.C.ex.max],
                                    eem.matrices[[i]][index.C.em.max-1, index.C.ex.max+1],
                                    eem.matrices[[i]][index.C.em.max-1, index.C.ex.max-1])

      if(C.intensity[i]==max.neighbor.intensity) {
        break
      } else {
        testable.values <- testable.values[!(testable.values %in% C.intensity[i])]
        testcount[i] <- testcount[i] + 1
      }
    }
    if(testcount[i]>20) { warning(paste0(sample.names[i], ": At least 20 failures to locate peak C maximum.")) }

    if(show.plots) {
      filled.contour(x=as.numeric(dimnames(eem.matrices[[i]])[[2]]), y=as.numeric(dimnames(eem.matrices[[i]])[[1]]),
                     z=t(eem.matrices[[i]]), color=matlab.like, nlevels=50,
                     plot.axes = {axis(1); axis(2); points(C.ex.max[i],C.em.max[i],pch=19)},
                     plot.title=title(main=eem.files[i]))
      notes[i] <- readline("Does this C peak maximum look OK? Type any notes, or just press ENTER.")
    }
  }

  indices <- data.frame(row.names=sample.names, BIX, FI, HIX, C.intensity, C.em.max, C.ex.max, ini.C.em.max, min.C.ex, testcount, notes)

  write.csv(indices, file=outputfile)

}
