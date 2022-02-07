library(ggplot2)
library(patchwork)
library(here)

# Get names of txt files
filenames <- dir(path = "example_traces/", pattern = "txt")

if (length(filenames)) {
  print(paste("Found", length(filenames), "traces."))
  } else {
  stop("No traces found")
  }

# Read in data
traces <- lapply(paste0("example_traces/", filenames), function(f) {
  t <- read.table(f, header = TRUE)
  names(t) <- c("Time", "Voltage")
  t
})

plot_all_traces <- function() {
  plot.list <- list()

  for (t in seq_len(traces)) {
    plot.list[[t]] <- ggplot(traces[[t]], aes(x = Time, y = Voltage)) +
      geom_line()
  }
  wrap_plots(plot.list, ncol = 3)
}

get_aps <- function(trace,
                    baseline_quantile = 0.7,
                    ap_quantile = 0.9,
                    min_ap_len = 30,
                    max_ap_len = 300,
                    verbose = FALSE) {
  # Get the action potentials positions in the trace
  # Parameters:
  #   - trace: the trace to get the action potentials positions from
  #   - baseline_quantile: the quantile used to define the baseline.
  #     Defaults to 0.7
  #   - ap_quantile: the quantile used to define the action potential.
  #     AP search begins when the trace crosses the baseline. For an
  #     event to be considered, the trace must cross this threshold as well.
  #     Defaults to 0.9
  #   - min_ap_len: the minimum length of an action potential, in ms.
  #     Defaults to 5
  #   - verbose: whether to print progress information. Defaults to FALSE

  # Determine baseline and AP threshold
  baseline <- quantile(trace$Voltage, baseline_quantile)
  ap_threshold <- quantile(trace$Voltage, ap_quantile)

  # Get the sampling interval
  sampling_int <- trace$Time[2] - trace$Time[1]
  # Convert min/max AP length to sampling points
  min_ap_len <- min_ap_len / sampling_int
  max_ap_len <- max_ap_len / sampling_int
  
  # Points where baseline is passed
  over_thr <- trace$Voltage > baseline
  # These are the points before and after the baseline
  changepoints <- which(diff(over_thr) != 0)

  # Remove initial changepoint, if trace starts above baseline
  if (trace$Voltage[1] > baseline) {
    changepoints <- changepoints[-1]
  }

  # Remove final changepoint, if trace ends above baseline
  if (trace$Voltage[length(trace$Voltage)] > baseline) {
    changepoints <- changepoints[1:(length(changepoints) - 1)]
  }

  changepoints <- matrix(changepoints, ncol = 2, byrow = TRUE)

  # Get the maximum between each pair of changepoints
  maxima <- apply(changepoints, 1, function(x) {
    which.max(trace$Voltage[x[1]:x[2]]) + x[1]
  })

  if (verbose) {
    cat(paste0("Found ", length(maxima), " events\n"))
  }

  # Get rid of all maxima that are below the AP threshold
  maxima <- maxima[trace$Voltage[maxima] > ap_threshold]

  if (verbose) {
    cat(paste0("Found ", length(maxima), " events above threshold\n"))
  }

  # Go back until we hit baseline to get the AP start
  start <- sapply(maxima, function(mx) {
    from <- max(1, mx - (max_ap_len %/% 2))
    
    res <- which.max(which(trace$Voltage[from:mx]<baseline)) + from
    
    if (!length(res))
      res <- mx - (max_ap_len %/% 2)
      
    res
  })

  # Go forward until we hit minimum to find the AHP
  stop <- apply(cbind(start, maxima), 1, function(sm) {
    to <- sm["maxima"] + max_ap_len - (sm["maxima"]-sm["start"])
  
    ahp_min <- which.min(trace$Voltage[sm["maxima"]:to]) + sm["maxima"]
    end <- which.max(which(trace$Voltage[ahp_min:to] < baseline)) + ahp_min

    if (!length(end))
      end <- to
    
    end ## We'll find the minimum later, after ensuring the max is ok
  })
  
  # We want to make sure APs don't overlap
  # Check if any start comes before the stop of the previous AP
  overlapping <- which(start[-1] < stop[1:length(stop)-1])
  stop[overlapping] <- start[overlapping+1] - 1
  
  minima <- apply(cbind(maxima, stop), 1, function(ms){
    res <- which.min(trace$Voltage[ms["maxima"]:ms["stop"]]) + ms["maxima"]
  })

  duration <- stop - start

  if (verbose) {
    cat(paste0("Event duration statistics",
                 "\n Min: ", min(duration) / sampling_int,
                 "\n Mean: ", mean(duration) / sampling_int,
                 "\n Max: ", max(duration) / sampling_int, "\n"))
  }

  to_keep <- which(duration > min_ap_len)

  if (verbose) {
    cat(paste("Removed", length(maxima) - length(to_keep), 
              "events below the minimum duration of", min_ap_len))
  }

  # Now get the APs
  aps <- lapply(seq_along(start), function(i){
    trace[start[i]:stop[i],]
  })
  
  res <- list(
    baseline = baseline,
    ap_threshold = ap_threshold,
    start = start[to_keep],
    stop = stop[to_keep],
    duration = duration[to_keep],
    min = minima[to_keep],
    max = maxima[to_keep],
    aps = aps
  )

  res
}

plot_aps <- function(trace, aps,
                     plot_max = TRUE, 
                     plot_min = TRUE,
                     plot_start = TRUE,
                     plot_stop = TRUE, 
                     plot_baseline = TRUE) {
  p <- ggplot(trace, aes(x = Time, y = Voltage)) +
    geom_line()

  if (plot_max) {
    p <- p + geom_point(
      data = trace[aps$max, ],
      color = "orange",
      cex = 3
    )
  }

    if (plot_min) {
    p <- p + geom_point(
      data = trace[aps$min, ],
      color = "blue",
      pch = 25,
      cex = 2
    )
  }
  
  if (plot_start) {
    p <- p + geom_point(
      data = trace[aps$start, ],
      color = "green",
      pch = 17,
      cex = 2.5
    )
  }

  if (plot_stop) {
    p <- p + geom_point(
      data = trace[aps$stop, ],
      color = "red",
      pch = 15,
      cex = 2.5
    )
  }

  if (plot_baseline) {
    p <- p +
      geom_hline(yintercept = aps$baseline, col = "red", lty = "dotted")
  }
  print(p)
  
  p
}

trace <- traces[[1]]

aps <- get_aps(trace,
  baseline_quantile = 0.5,
  max_ap_len = 500,
  verbose = TRUE
)

plot_aps(trace, aps) 

p <- lapply(aps$aps, function(ap){
  ggplot(ap, aes(x = Time-min(Time), y = Voltage)) + geom_line()
  })

wrap_plots(p)


##### PCA ANALYSIS #####

# Now we want to align spikes on the maximum
dist_to_max <- sapply(aps$aps, function(x) {
  max_pos <- which.max(x$Voltage)
  max_to_end <- nrow(x) - max_pos
  
  c(start_to_max = max_pos, max_to_end = max_to_end)
})

max_left <- max(dist_to_max[1,])
max_right <- max(dist_to_max[2,])

padded_aps <- sapply(1:length(aps$aps), function(i) {
  ap <- aps$aps[[i]]$Voltage
  ap <- c(rep(ap[1], max_left - dist_to_max[1,i]), ap)
  ap <- c(ap, rep(ap[length(ap)], max_right - dist_to_max[2,i]))
  
  ap
})

padded_aps <- data.frame(padded_aps)

p <- apply(padded_aps, 2, function(ap){
  ggplot(data.frame(Sample = seq_along(ap), Voltage = ap), aes(x = Sample, y = Voltage)) + 
    geom_line()
})

wrap_plots(p)

