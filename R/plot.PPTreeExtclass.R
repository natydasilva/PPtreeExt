#' Plot Projection Pursuit Classification Tree
#'
#' @description
#' Visualizes a Projection Pursuit (PP) classification tree  using grid graphics.
#' The function creates a hierarchical tree diagram showing the structure of 
#' splits and terminal nodes with class assignments. Supports automatic scaling
#' for large trees.
#'
#' @param x An object of class \code{PPtreeclass} or \code{PPtreeExtclass} 
#'   containing the tree structure, projection vectors, and split points.
#' @param font.size Numeric. Font size for text labels in the plot. 
#'   Default is 17. Will be automatically reduced for large trees when 
#'   \code{auto.scale = TRUE}.
#' @param width.size Numeric. Width scaling factor for graphical elements 
#'   (nodes, edges). Default is 1. Will be automatically adjusted for large 
#'   trees when \code{auto.scale = TRUE}.
#' @param main Character string. Main title for the plot. 
#'   Default is "Projection Pursuit Classification Tree".
#' @param sub Character string or NULL. Subtitle for the plot. 
#'   Default is NULL (no subtitle).
#' @param auto.scale Logical. If TRUE (default), automatically adjusts plot 
#'   dimensions and font size based on tree size. Recommended for trees with 
#'   more than 20 terminal nodes or depth greater than 15.
#' @param min.width Numeric or NULL. Minimum width (number of terminal nodes) 
#'   for the plot when \code{auto.scale = FALSE}. If NULL, uses the number of 
#'   classes in the data. Default is NULL.
#' @param min.height Numeric or NULL. Minimum height (tree depth) for the plot 
#'   when \code{auto.scale = FALSE}. If NULL, uses calculated tree depth. 
#'   Default is NULL.
#' @param ... Additional arguments (currently not used).
#'
#' @details
#' The plot displays:
#' \itemize{
#'   \item \strong{Internal nodes}: Shown as ellipses with the projection used 
#'     for splitting (e.g., "proj1 * X")
#'   \item \strong{Terminal nodes}: Shown as gray rectangles with the assigned 
#'     class label
#'   \item \strong{Edges}: Labeled with split rules ("< cutN" for left child, 
#'     ">= cutN" for right child)
#'   \item \strong{Node IDs}: Small boxes at the top of each node
#' }
#'
#' \strong{Auto-scaling behavior:}
#' When \code{auto.scale = TRUE} and the tree has more than 20 terminal nodes 
#' or depth greater than 15:
#' \itemize{
#'   \item Font size is reduced: \code{max(10, 17 - floor((n_terminal - 20) / 5))}
#'   \item Width size is reduced: \code{max(0.7, 1 - (n_terminal - 20) * 0.01)}
#'   \item Width is set to: \code{max(n_terminal, n_classes)}
#'   \item Height is set to: tree depth
#' }
#'
#' \strong{Manual scaling:}
#' When \code{auto.scale = FALSE}, you can control dimensions with 
#' \code{min.width} and \code{min.height}.
#'
#' @return 
#' Invisibly returns a list with:
#' \item{width}{Numeric. The width used for plotting (number of terminal nodes)}
#' \item{height}{Numeric. The height used for plotting (tree depth)}
#' \item{font.size}{Numeric. The final font size used (after auto-scaling)}
#'
#' @note
#' This function requires the \code{grid} package. It will create a new graphics 
#' page using \code{grid.newpage()}.
#'
#' For very large trees (>50 terminal nodes), consider:
#' \itemize{
#'   \item Exporting to a large PNG or PDF file
#'   \item Manually reducing \code{font.size} further
#'   \item Pruning the tree before plotting
#' }

#'
#' @examples
#' \dontrun{
#' library(grid)
#' # Example with penguins dataset
#' data(penguins)
#' penguins <- na.omit(penguins[, -c(2,7)])
#' penguins_ppt <- PPtreeExtclass(species~bill_len + bill_dep +flipper_len +
#'  body_mass,  data = penguins, PPmethod = "PDA", srule = FALSE )
#' 
#' plot(penguins_ppt,
#'      main = "Penguins Classification with PPtreeExt",
#'       font.size = 8,        
#' width.size = 0.7)
#' }
#'
#' @export
#' @import grid
#' @method plot PPtreeExtclass
#' 
plot.PPtreeExtclass <- function(x, font.size = 17, width.size = 1, 
                             main = "Projection Pursuit Classification Tree", 
                             sub = NULL, auto.scale = TRUE, 
                             min.width = NULL, min.height = NULL, ...) {
  
  PPtreeobj <- x
  
  # ===== auxiliary =====
  
  # compute tree depth
  
  calc.depth <- function(PPtreeobj) {
    TS <- PPtreeobj$Tree.Struct
    i <- 1
    flag.L <- rep(FALSE, nrow(TS))
    keep.track <- 1
    depth.track <- 0
    depth <- 0
    
    while (sum(flag.L) != nrow(TS)) {
      if (!flag.L[i]) {
        if (TS[i, 2] == 0) {
          flag.L[i] <- TRUE
          id.l <- length(keep.track) - 1
          i <- keep.track[id.l]
          depth <- depth - 1
        } else if (!flag.L[TS[i, 2]]) {
          depth <- depth + 1
          i <- TS[TS[i, 2], 1]
        } else {
          depth <- depth + 1
          flag.L[i] <- TRUE
          i <- TS[TS[i, 3], 1]
        }
        keep.track <- c(keep.track, i)
        depth.track <- c(depth.track, depth)
      } else {
        id.l <- id.l - 1
        i <- keep.track[id.l]
        depth <- depth.track[id.l]
      }
    }
    depth <- max(depth.track) + 2
    return(depth)
  }
  
  # Count final nodes in specific directions 
  n.final <- function(PPtreeobj, node.id, direction) {
    TS <- PPtreeobj$Tree.Struct
    n.leaf <- 0
    
    if (direction == "left") {
      keep.id <- TS[node.id, 2]
    } else if (direction == "right") {
      keep.id <- TS[node.id, 3]
    }
    
    i <- 1
    while (i <= length(keep.id)) {
      if (TS[keep.id[i], 2] == 0) {
        n.leaf <- n.leaf + 1
        i <- i + 1
      } else {
        keep.id <- c(keep.id, TS[keep.id[i], 2:3])
        i <- i + 1
      }
    }
    return(n.leaf)
  }
  
  # Count total final noes
  count.terminal.nodes <- function(PPtreeobj) {
    TS <- PPtreeobj$Tree.Struct
    sum(TS[, 2] == 0)
  }
  
  # labels of edges
  edge.lable.PPtree <- function(PPtreeobj, node.id, left = TRUE) {
    TS <- PPtreeobj$Tree.Struct
    if (left) {
      text.t <- paste("< cut", TS[node.id, 4], sep = "")
    } else {
      text.t <- paste(">= cut", TS[node.id, 4], sep = "")
    }
    
    grid.rect(gp = gpar(fill = "white", lty = 1, col = "grey95"), 
              width = unit(width.size, "strwidth", text.t) * 1.2)
    grid.text(text.t, just = "center", gp = gpar(fontsize = font.size))
  }
  
  # internal nodes
  node.inner.PPtree <- function(PPtreeobj, node.id) {
    TS <- PPtreeobj$Tree.Struct
    
    label.t <- paste("proj", TS[node.id, 4], " * X", sep = "")
    
    Inner.Node.V <- viewport(
      x = unit(0.5, "npc"), 
      y = unit(0.5, "npc"), 
      width = unit(width.size * 1.5, "strwidth", label.t), 
      height = unit(width.size * 2, "lines")
    )
    pushViewport(Inner.Node.V)
    
    # ellipse drawing
    xell <- c(seq(0, 0.2, by = 0.01), seq(0.2, 0.8, by = 0.05), seq(0.8, 1, by = 0.01))
    yell <- sqrt(xell * (1 - xell))
    grid.polygon(
      x = unit(c(xell, rev(xell)), "npc"), 
      y = unit(c(yell, -yell) + 0.5, "npc"), 
      gp = gpar(fill = "white")
    )
    
    grid.text(label.t, y = 0.3, gp = gpar(fontsize = font.size), 
              just = c("center", "bottom"))
    
    # node ID
    Inner.Node.Id.V <- viewport(
      x = unit(0.5, "npc"), 
      y = unit(1, "npc"), 
      width = max(unit(1, "lines"), unit(1.2, "strwidth", as.character(node.id))), 
      height = max(unit(1, "lines"), unit(1.2, "strheight", as.character(node.id))), 
      just = c("center", "center"), 
      gp = gpar(fontsize = font.size)
    )
    pushViewport(Inner.Node.Id.V)
    grid.rect(gp = gpar(fill = "white", lty = "solid", fontsize = font.size))
    grid.text(node.id, gp = gpar(fontsize = font.size))
    popViewport()
    upViewport()
  }
  
  # Nodo terminal
  node.terminal.PPtree <- function(PPtreeobj, node.id) {
    TS <- PPtreeobj$Tree.Struct
    gName <- names(table(PPtreeobj$origclass))
    
    # Verify that class index is valid
    class.idx <- TS[node.id, 3]
    if (is.na(class.idx) || class.idx < 1 || class.idx > length(gName)) {
      gN <- "NA"
    } else {
      gN <- gName[class.idx]
    }
    
    temp <- strsplit(as.character(gN), split = "")[[1]]
    gN.width <- length(temp)
    set.unit <- (sum(tolower(temp) != temp) * 0.65 + 
                   sum(tolower(temp) == temp) * 0.5) / gN.width
    
    Terminal.Node.V <- viewport(
      x = unit(0.5, "npc"), 
      y = unit(0.8, "npc"), 
      width = unit(1, "lines") * 1.5, 
      height = unit(set.unit, "lines") * (gN.width + 2), 
      just = c("center", "top")
    )
    pushViewport(Terminal.Node.V)
    
    grid.rect(gp = gpar(fill = "lightgray"))
    grid.text(y = 0.05, gN, gp = gpar(fontsize = font.size), 
              rot = 90, just = "left")
    
    # ID del nodo
    Terminal.Node.Id.V <- viewport(
      x = unit(0.5, "npc"), 
      y = unit(1, "npc"), 
      width = max(unit(1, "lines"), unit(1.2, "strwidth", as.character(node.id))), 
      height = max(unit(1, "lines"), unit(1.2, "strheight", as.character(node.id))), 
      just = c("center", "center"), 
      gp = gpar(fontsize = font.size)
    )
    pushViewport(Terminal.Node.Id.V)
    grid.rect(gp = gpar(fill = "lightgray", lty = "solid", fontsize = font.size))
    grid.text(node.id, gp = gpar(fontsize = font.size))
    popViewport()
    upViewport()
  }
  
  # Draw the tree recursively
  plotPPtree <- function(PPtreeobj, node.id, xlim, ylim) {
    TS <- PPtreeobj$Tree.Struct
    
    # terminal node
    if (TS[node.id, 2] == 0) {
      x <- xlim[1] + 0.5
      y <- ylim[2] - 1
      
      Final.Node.V <- viewport(
        x = unit(x, "native"), 
        y = unit(y, "native"), 
        width = unit(1, "native"), 
        height = unit(1, "native") - unit(2, "lines"), 
        just = c("center", "top")
      )
      pushViewport(Final.Node.V)
      node.terminal.PPtree(PPtreeobj, node.id)
      upViewport()
      return(NULL)
    }
    
    # Internal node
    nl <- n.final(PPtreeobj, node.id, "left")
    nr <- n.final(PPtreeobj, node.id, "right")
    x0 <- xlim[1] + nl
    y0 <- max(ylim) - 1
    
    lf <- ifelse(TS[TS[node.id, 2], 2] == 0, 0.5, 
                 n.final(PPtreeobj, TS[node.id, 2], "right"))
    rf <- ifelse(TS[TS[node.id, 3], 2] == 0, 0.5, 
                 n.final(PPtreeobj, TS[node.id, 3], "left"))
    
    x1l <- x0 - lf
    x1r <- x0 + rf
    y1 <- y0 - 1
    
    # Draw lines
    grid.lines(x = unit(c(x0, x1l), "native"), y = unit(c(y0, y1), "native"))
    grid.lines(x = unit(c(x0, x1r), "native"), y = unit(c(y0, y1), "native"))
    
    # Draw nodes
    node.V <- viewport(
      x = unit(x0, "native"), 
      y = unit(y0, "native"), 
      width = unit(1, "native"), 
      height = unit(1, "native") - unit(1, "lines")
    )
    pushViewport(node.V)
    node.inner.PPtree(PPtreeobj, node.id)
    upViewport()
    
    # Label edges
    ylpos <- y0 - 0.6
    yrpos <- y0 - 0.45
    xlpos <- x0 - (x0 - x1l) * 0.6
    xrpos <- x0 - (x0 - x1r) * 0.45
    
    LeftEdge.V <- viewport(
      x = unit(xlpos, "native"), 
      y = unit(ylpos, "native"), 
      width = unit(abs(xlpos - xrpos), "native"), 
      height = unit(1, "lines") * 1.2
    )
    pushViewport(LeftEdge.V)
    edge.lable.PPtree(PPtreeobj, node.id, left = TRUE)
    upViewport()
    
    RightEdge.V <- viewport(
      x = unit(xrpos, "native"), 
      y = unit(yrpos, "native"), 
      width = unit(abs(xlpos - xrpos), "native"), 
      height = unit(1, "lines")
    )
    pushViewport(RightEdge.V)
    edge.lable.PPtree(PPtreeobj, node.id, left = FALSE)
    upViewport()
    
    # recursion for kids
    plotPPtree(PPtreeobj, TS[node.id, 2], c(xlim[1], x0), c(1, y1 + 1))
    plotPPtree(PPtreeobj, TS[node.id, 3], c(x0, xlim[2]), c(1, y1 + 1))
  }
  
  # ===== Compute dimensions =====
  
  # number of terminal nodes (width of the tree)
  n_terminal <- count.terminal.nodes(PPtreeobj)
  
  # Depth of the tree  (height)
  depth <- calc.depth(PPtreeobj)
  
  # Automatic dimension scaling
  if (auto.scale) {
    # For big trees adjust width 
    nx <- max(n_terminal, length(table(PPtreeobj$origclass)))
    
    # Fot depth trees adjust height 
    ny <- depth
    
    # Adjust fornt size if the tree is  big 
    if (n_terminal > 20 || depth > 15) {
      font.size <- max(10, 17 - floor((n_terminal - 20) / 5))
      width.size <- max(0.7, 1 - (n_terminal - 20) * 0.01)
      cat("Auto-scaling: Terminal nodes =", n_terminal, 
          ", Depth =", depth, 
          ", Font size =", font.size, "\n")
    }
  } else {
    nx <- ifelse(is.null(min.width), length(table(PPtreeobj$origclass)), min.width)
    ny <- ifelse(is.null(min.height), depth, min.height)
  }
  
  # Asure min dimensions
  nx <- max(nx, 3)
  ny <- max(ny, 3)
  
  cat("Tree dimensions: Width =", nx, ", Height =", ny, "\n")
  
  # ===== Plot tree =====
  
  grid.newpage()
  
  PPtree.Main.V <- viewport(
    layout = grid.layout(
      3, 3, 
      heights = unit(c(3, 1, 1), c("lines", "null", "lines")), 
      widths = unit(c(1, 1, 1), c("lines", "null", "lines"))
    )
  )
  pushViewport(PPtree.Main.V)
  
  # Title
  PPtree.title.V <- viewport(layout.pos.col = 2, layout.pos.row = 1)
  pushViewport(PPtree.title.V)
  grid.text(
    y = unit(1, "lines"), 
    paste("\n", main, sep = ""), 
    just = "center", 
    gp = gpar(fontsize = font.size)
  )
  upViewport()
  
  # Tree
  PPtree.Tree.V <- viewport(
    layout.pos.col = 2, 
    layout.pos.row = 2, 
    xscale = c(0, nx), 
    yscale = c(0, ny + 1)
  )
  pushViewport(PPtree.Tree.V)
  
  plotPPtree(PPtreeobj, 1, c(0, nx), ylim = c(1, ny + 1))
  
  # Subtitle
  if (!is.null(sub)) {
    grid.text(
      y = unit(1, "lines"), 
      sub, 
      just = "center", 
      gp = gpar(fontsize = font.size * 0.7)
    )
  }
  
  popViewport(2)
  
  invisible(list(width = nx, height = ny, font.size = font.size))
}
plot.PPtreeclass <- function(x, font.size = 17, width.size = 1, 
                             main = "Projection Pursuit Classification Tree", 
                             sub = NULL, auto.scale = TRUE, 
                             min.width = NULL, min.height = NULL, ...) {
    
    PPtreeobj <- x
    
    # ===== Auxiliary functions =====
    
    # Compute the depth of the tree
    calc.depth <- function(PPtreeobj) {
        TS <- PPtreeobj$Tree.Struct
        i <- 1
        flag.L <- rep(FALSE, nrow(TS))
        keep.track <- 1
        depth.track <- 0
        depth <- 0
        
        while (sum(flag.L) != nrow(TS)) {
            if (!flag.L[i]) {
                if (TS[i, 2] == 0) {
                    flag.L[i] <- TRUE
                    id.l <- length(keep.track) - 1
                    i <- keep.track[id.l]
                    depth <- depth - 1
                } else if (!flag.L[TS[i, 2]]) {
                    depth <- depth + 1
                    i <- TS[TS[i, 2], 1]
                } else {
                    depth <- depth + 1
                    flag.L[i] <- TRUE
                    i <- TS[TS[i, 3], 1]
                }
                keep.track <- c(keep.track, i)
                depth.track <- c(depth.track, depth)
            } else {
                id.l <- id.l - 1
                i <- keep.track[id.l]
                depth <- depth.track[id.l]
            }
        }
        depth <- max(depth.track) + 2
        return(depth)
    }
    
    # Count the final nodes in specific directions
    n.final <- function(PPtreeobj, node.id, direction) {
        TS <- PPtreeobj$Tree.Struct
        n.leaf <- 0
        
        if (direction == "left") {
            keep.id <- TS[node.id, 2]
        } else if (direction == "right") {
            keep.id <- TS[node.id, 3]
        }
        
        i <- 1
        while (i <= length(keep.id)) {
            if (TS[keep.id[i], 2] == 0) {
                n.leaf <- n.leaf + 1
                i <- i + 1
            } else {
                keep.id <- c(keep.id, TS[keep.id[i], 2:3])
                i <- i + 1
            }
        }
        return(n.leaf)
    }
    
    # Count the total terminal nodes
    count.terminal.nodes <- function(PPtreeobj) {
        TS <- PPtreeobj$Tree.Struct
        sum(TS[, 2] == 0)
    }
    
    # Label of edges
    edge.lable.PPtree <- function(PPtreeobj, node.id, left = TRUE) {
        TS <- PPtreeobj$Tree.Struct
        if (left) {
            text.t <- paste("< cut", TS[node.id, 4], sep = "")
        } else {
            text.t <- paste(">= cut", TS[node.id, 4], sep = "")
        }
        
        grid.rect(gp = gpar(fill = "white", lty = 1, col = "grey95"), 
                  width = unit(width.size, "strwidth", text.t) * 1.2)
        grid.text(text.t, just = "center", gp = gpar(fontsize = font.size))
    }
    
    # Internal node
    node.inner.PPtree <- function(PPtreeobj, node.id) {
        TS <- PPtreeobj$Tree.Struct
        
        label.t <- paste("proj", TS[node.id, 4], " * X", sep = "")
        
        Inner.Node.V <- viewport(
            x = unit(0.5, "npc"), 
            y = unit(0.5, "npc"), 
            width = unit(width.size * 1.5, "strwidth", label.t), 
            height = unit(width.size * 2, "lines")
        )
        pushViewport(Inner.Node.V)
        
        # Plot elliose
        xell <- c(seq(0, 0.2, by = 0.01), seq(0.2, 0.8, by = 0.05), seq(0.8, 1, by = 0.01))
        yell <- sqrt(xell * (1 - xell))
        grid.polygon(
            x = unit(c(xell, rev(xell)), "npc"), 
            y = unit(c(yell, -yell) + 0.5, "npc"), 
            gp = gpar(fill = "white")
        )
        
        grid.text(label.t, y = 0.3, gp = gpar(fontsize = font.size), 
                  just = c("center", "bottom"))
        
        # node ID
        Inner.Node.Id.V <- viewport(
            x = unit(0.5, "npc"), 
            y = unit(1, "npc"), 
            width = max(unit(1, "lines"), unit(1.2, "strwidth", as.character(node.id))), 
            height = max(unit(1, "lines"), unit(1.2, "strheight", as.character(node.id))), 
            just = c("center", "center"), 
            gp = gpar(fontsize = font.size)
        )
        pushViewport(Inner.Node.Id.V)
        grid.rect(gp = gpar(fill = "white", lty = "solid", fontsize = font.size))
        grid.text(node.id, gp = gpar(fontsize = font.size))
        popViewport()
        upViewport()
    }
    
    # Terminal node
    node.terminal.PPtree <- function(PPtreeobj, node.id) {
        TS <- PPtreeobj$Tree.Struct
        gName <- names(table(PPtreeobj$origclass))
        
        # Verify that class index is valid
        class.idx <- TS[node.id, 3]
        if (is.na(class.idx) || class.idx < 1 || class.idx > length(gName)) {
            gN <- "NA"
        } else {
            gN <- gName[class.idx]
        }
        
        temp <- strsplit(as.character(gN), split = "")[[1]]
        gN.width <- length(temp)
        set.unit <- (sum(tolower(temp) != temp) * 0.65 + 
                     sum(tolower(temp) == temp) * 0.5) / gN.width
        
        Terminal.Node.V <- viewport(
            x = unit(0.5, "npc"), 
            y = unit(0.8, "npc"), 
            width = unit(1, "lines") * 1.5, 
            height = unit(set.unit, "lines") * (gN.width + 2), 
            just = c("center", "top")
        )
        pushViewport(Terminal.Node.V)
        
        grid.rect(gp = gpar(fill = "lightgray"))
        grid.text(y = 0.05, gN, gp = gpar(fontsize = font.size), 
                  rot = 90, just = "left")
        
        # Node ID
        Terminal.Node.Id.V <- viewport(
            x = unit(0.5, "npc"), 
            y = unit(1, "npc"), 
            width = max(unit(1, "lines"), unit(1.2, "strwidth", as.character(node.id))), 
            height = max(unit(1, "lines"), unit(1.2, "strheight", as.character(node.id))), 
            just = c("center", "center"), 
            gp = gpar(fontsize = font.size)
        )
        pushViewport(Terminal.Node.Id.V)
        grid.rect(gp = gpar(fill = "lightgray", lty = "solid", fontsize = font.size))
        grid.text(node.id, gp = gpar(fontsize = font.size))
        popViewport()
        upViewport()
    }
    
    # Plot the tree recursively
    plotPPtree <- function(PPtreeobj, node.id, xlim, ylim) {
        TS <- PPtreeobj$Tree.Struct
        
        # Terminal node
        if (TS[node.id, 2] == 0) {
            x <- xlim[1] + 0.5
            y <- ylim[2] - 1
            
            Final.Node.V <- viewport(
                x = unit(x, "native"), 
                y = unit(y, "native"), 
                width = unit(1, "native"), 
                height = unit(1, "native") - unit(2, "lines"), 
                just = c("center", "top")
            )
            pushViewport(Final.Node.V)
            node.terminal.PPtree(PPtreeobj, node.id)
            upViewport()
            return(NULL)
        }
        
        # Internal node
        nl <- n.final(PPtreeobj, node.id, "left")
        nr <- n.final(PPtreeobj, node.id, "right")
        x0 <- xlim[1] + nl
        y0 <- max(ylim) - 1
        
        lf <- ifelse(TS[TS[node.id, 2], 2] == 0, 0.5, 
                     n.final(PPtreeobj, TS[node.id, 2], "right"))
        rf <- ifelse(TS[TS[node.id, 3], 2] == 0, 0.5, 
                     n.final(PPtreeobj, TS[node.id, 3], "left"))
        
        x1l <- x0 - lf
        x1r <- x0 + rf
        y1 <- y0 - 1
        
        # Dibuja líneas
        grid.lines(x = unit(c(x0, x1l), "native"), y = unit(c(y0, y1), "native"))
        grid.lines(x = unit(c(x0, x1r), "native"), y = unit(c(y0, y1), "native"))
        
        # Dibuja nodo
        node.V <- viewport(
            x = unit(x0, "native"), 
            y = unit(y0, "native"), 
            width = unit(1, "native"), 
            height = unit(1, "native") - unit(1, "lines")
        )
        pushViewport(node.V)
        node.inner.PPtree(PPtreeobj, node.id)
        upViewport()
        
        # Label edges
        ylpos <- y0 - 0.6
        yrpos <- y0 - 0.45
        xlpos <- x0 - (x0 - x1l) * 0.6
        xrpos <- x0 - (x0 - x1r) * 0.45
        
        LeftEdge.V <- viewport(
            x = unit(xlpos, "native"), 
            y = unit(ylpos, "native"), 
            width = unit(abs(xlpos - xrpos), "native"), 
            height = unit(1, "lines") * 1.2
        )
        pushViewport(LeftEdge.V)
        edge.lable.PPtree(PPtreeobj, node.id, left = TRUE)
        upViewport()
        
        RightEdge.V <- viewport(
            x = unit(xrpos, "native"), 
            y = unit(yrpos, "native"), 
            width = unit(abs(xlpos - xrpos), "native"), 
            height = unit(1, "lines")
        )
        pushViewport(RightEdge.V)
        edge.lable.PPtree(PPtreeobj, node.id, left = FALSE)
        upViewport()
        
        # Child recursion
        plotPPtree(PPtreeobj, TS[node.id, 2], c(xlim[1], x0), c(1, y1 + 1))
        plotPPtree(PPtreeobj, TS[node.id, 3], c(x0, xlim[2]), c(1, y1 + 1))
    }
    
    # ===== COMPUTE DIMENSIONS =====
    
    # Number of terminal nodes (width of the tree)
    n_terminal <- count.terminal.nodes(PPtreeobj)
    
    # Depth of the tree (height)
    depth <- calc.depth(PPtreeobj)
    
    # Automatic adjust of dimensions
    if (auto.scale) {
        # For big trees, adjust the width
        nx <- max(n_terminal, length(table(PPtreeobj$origclass)))
        
        # For deepth trees, adjust height 
        ny <- depth
        
        # Adjust font size if the tree is big
        if (n_terminal > 20 || depth > 15) {
            font.size <- max(10, 17 - floor((n_terminal - 20) / 5))
            width.size <- max(0.7, 1 - (n_terminal - 20) * 0.01)
            cat("Auto-scaling: Terminal nodes =", n_terminal, 
                ", Depth =", depth, 
                ", Font size =", font.size, "\n")
        }
    } else {
        nx <- ifelse(is.null(min.width), length(table(PPtreeobj$origclass)), min.width)
        ny <- ifelse(is.null(min.height), depth, min.height)
    }
    
    # Minimal dimensions
    nx <- max(nx, 3)
    ny <- max(ny, 3)
    
    cat("Tree dimensions: Width =", nx, ", Height =", ny, "\n")
    
    # ===== Drow the treeL =====
    
    grid.newpage()
    
    PPtree.Main.V <- viewport(
        layout = grid.layout(
            3, 3, 
            heights = unit(c(3, 1, 1), c("lines", "null", "lines")), 
            widths = unit(c(1, 1, 1), c("lines", "null", "lines"))
        )
    )
    pushViewport(PPtree.Main.V)
    
    # Title
    PPtree.title.V <- viewport(layout.pos.col = 2, layout.pos.row = 1)
    pushViewport(PPtree.title.V)
    grid.text(
        y = unit(1, "lines"), 
        paste("\n", main, sep = ""), 
        just = "center", 
        gp = gpar(fontsize = font.size)
    )
    upViewport()
    
    # Tree
    PPtree.Tree.V <- viewport(
        layout.pos.col = 2, 
        layout.pos.row = 2, 
        xscale = c(0, nx), 
        yscale = c(0, ny + 1)
    )
    pushViewport(PPtree.Tree.V)
    
    plotPPtree(PPtreeobj, 1, c(0, nx), ylim = c(1, ny + 1))
    
    # Subtitle
    if (!is.null(sub)) {
        grid.text(
            y = unit(1, "lines"), 
            sub, 
            just = "center", 
            gp = gpar(fontsize = font.size * 0.7)
        )
    }
    
    popViewport(2)
    
    invisible(list(width = nx, height = ny, font.size = font.size))
}