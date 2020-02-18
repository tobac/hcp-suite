# standard gifti image plotting functions, intended for calling when plotting in knitr files (or similar). 
# see http://mvpa.blogspot.com/2018/06/tutorial-plotting-gifti-images-in-r.html for more explanation
# code written by Joset A. Etzel (jetzel@wustl.edu) and may be adapted, provided this source is cited. 1 June 2018.
# Cognitive Control & Psychopathology Lab, Psychological & Brain Sciences, Washington University in St. Louis (USA)
# last updated 10 April 2019 by Jo Etzel.

# colors to use for positive (warm) and negative (cool) values, generated using http://www.perbang.dk/rgbgradient/
cols.warm <- c("#E50800", "#E50C00", "#E51001", "#E61402", "#E61803", "#E61C03", "#E72004", "#E72405", "#E72806", "#E82D06", 
               "#E83107", "#E83508", "#E93909", "#E93D09", "#E9410A", "#EA450B", "#EA490C", "#EB4E0C", "#EB520D", "#EB560E", 
               "#EC5A0F", "#EC5E10", "#EC6210", "#ED6611", "#ED6A12", "#ED6E13", "#EE7313", "#EE7714", "#EE7B15", "#EF7F16", 
               "#EF8316", "#F08717", "#F08B18", "#F08F19", "#F19419", "#F1981A", "#F19C1B", "#F2A01C", "#F2A41C", "#F2A81D", 
               "#F3AC1E", "#F3B01F", "#F3B420", "#F4B920", "#F4BD21", "#F5C122", "#F5C523", "#F5C923", "#F6CD24", "#F6D125", 
               "#F6D526", "#F7DA26", "#F7DE27", "#F7E228", "#F8E629", "#F8EA29", "#F8EE2A", "#F9F22B", "#F9F62C", "#FAFA2D");

cols.cool <- c("#001CE5", "#001FE4", "#0123E4", "#0227E4", "#032BE3", "#032EE3", "#0432E3", "#0536E3", "#063AE2", "#063EE2", 
               "#0741E2", "#0845E2", "#0949E1", "#094DE1", "#0A50E1", "#0B54E1", "#0C58E0", "#0C5CE0", "#0D60E0", "#0E63E0", 
               "#0F67DF", "#106BDF", "#106FDF", "#1172DF", "#1276DE", "#137ADE", "#137EDE", "#1482DE", "#1585DD", "#1689DD", 
               "#168DDD", "#1791DD", "#1894DC", "#1998DC", "#199CDC", "#1AA0DC", "#1BA4DB", "#1CA7DB", "#1CABDB", "#1DAFDB", 
               "#1EB3DA", "#1FB6DA", "#20BADA", "#20BEDA", "#21C2D9", "#22C6D9", "#23C9D9", "#23CDD9", "#24D1D8", "#25D5D8", 
               "#26D8D8", "#26DCD8", "#27E0D7", "#28E4D7", "#29E8D7", "#29EBD7", "#2AEFD6", "#2BF3D6", "#2CF7D6", "#2DFBD6");


# adapted from gifti package gifti_map_val function. See notes in "D:\gitFiles_ccplabwustl\R01\Jo\for800msecTR\gifti_test.R"
gifti.map <- function (pointset, triangle, values) {  # pointset <- surf.L$data$pointset; triangle <- surf.L$data$triangle; values <- over.L$data[[1]][,1]; 
  tri <- as.vector(t(triangle) + 1);  # transpose, unwrap + 1; +1 for 0-based (input) vs. 1-based (R)  
  out.vals <- values[tri[seq(from=1, to=length(tri), by=3)]];  # color triangle by FIRST point in each:
  pointset <- pointset[tri, ];
  
  return(list(coords=pointset, values=out.vals))
}


# surface plotting function
# gii is the object returned by gifti.map
# ttl is a string to display at the plot title (use "" for no title)
# pos.lims and neg.lims are the limits for the positive and negative (respectively) color scaling. Values from 0 to pos.lims[1] are not plotted,
# nor are values from neg.lims[2] to 0. Values larger than pos.lims[2] or smaller than neg.lims[1] are shown in the hottest or coolest color.
# print.scaling.label controls whether the color range label is written on the finished plot.
# which.surface indicates which underlay surface anatomy is being plotted, HCP or fsaverage5.

plot_surface <- function(gii, ttl, pos.lims=c(NA,NA), neg.lims=c(NA,NA), print.scaling.label=TRUE, which.surface="HCP") {  
  # gii <- tmp; ttl <- paste(lbls[i], "L", sess.ids[ssid]); pos.lims <- c(2,6); neg.lims <- c(-6,-2);
  full.min <- min(gii$values, na.rm=TRUE);  # store original minimum and maximum for later text label
  full.max <- max(gii$values, na.rm=TRUE);
  sc.text <- "";   # start scaling label with a blank string
  
  # first fix the values for proper plotting
  if (!is.na(pos.lims[1])) {    # showing positive values, so fix too-big values (so don't disappear when plotting)
    inds <- which(gii$values > pos.lims[2]);   
    if (length(inds) > 0) { gii$values[inds] <- pos.lims[2]; }
  }
  if (!is.na(neg.lims[1])) {    # showing negative values, so fix too-small values
    inds <- which(gii$values < neg.lims[1]);   
    if (length(inds) > 0) { gii$values[inds] <- neg.lims[1]; }
  }
  if (!is.na(pos.lims[1]) & !is.na(neg.lims[1])) {    # both positive and negtive, so fix middle (zero-ish) values
    inds <- which(gii$values > neg.lims[2] & gii$values < pos.lims[1]);   # NA out middle values so don't get plotted
    if (length(inds) > 0) { gii$values[inds] <- NA; }
    c.lims <- c(neg.lims[1], pos.lims[2]);   # set color scaling to smallest negative to biggest positive
    cols <- c(rev(cols.cool), cols.warm);   # both hot and cool colors
    sc.text <- paste("colors:", neg.lims[1], "to", neg.lims[2], "&", pos.lims[1], "to", pos.lims[2]);
  }
  if (!is.na(pos.lims[1]) & is.na(neg.lims[1])) {   # only positive values
    c.lims <- pos.lims;
    cols <- cols.warm;       
    sc.text <- paste("colors:", pos.lims[1], "to", pos.lims[2]);
  }
  if (is.na(pos.lims[1]) & !is.na(neg.lims[1])) {   # only negative values
    c.lims <- neg.lims;
    cols <- rev(cols.cool);
    sc.text <- paste("colors:", neg.lims[1], "to", neg.lims[2]);
  }
  
  # slightly different z-axis scaling and title positioning seems to work better for HCP and fsaverage5 (fMRIPrep) surface anatomies
  if (which.surface == "HCP") { z.lim <- c(-60,90); ttl.line <- -2.5; }
  if (which.surface == "fsaverage5") { z.lim <- c(-130,90); ttl.line <- -1; }
  
  # actually plot
  triangle3D(tri=gii$coords, theta=90, phi=0, ltheta=90, lphi=0, bty='n', colkey=FALSE, zlim=z.lim, d=6, colvar=gii$values, col=cols, clim=c.lims, facets=TRUE, resfac=0.01);  
  mtext(text=ttl, side=3, line=ttl.line, cex=0.6, adj=0);
  
  
  triangle3D(tri=gii$coords, theta=270, phi=0, ltheta=270, lphi=0, bty='n', colkey=FALSE, zlim=z.lim, d=6, colvar=gii$values, col=cols, clim=c.lims, facets=TRUE, resfac=0.01);  
  if (print.scaling.label == TRUE) {
    mtext(text=paste0("range: ", round(full.min,3), " to ", round(full.max,3), "\n", sc.text), side=1, line=-2.5, cex=0.5, adj=1);
  }
  if (print.scaling.label == "flipped") {
    mtext(text=paste0(sc.text, "\nrange: ", round(full.min,3), " to ", round(full.max,3)), side=1, line=-2.5, cex=0.5, adj=1);
  }
}
