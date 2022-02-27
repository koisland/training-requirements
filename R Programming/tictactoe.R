library(R6)


TicTacToe <- R6Class("TicTacToe",
  public = list(
    initialize = function(empty_char, chars, grid_size) {
      
      # Initialize game.
      # @param: empty_char - character to use as empty space on grid
      # @param: chars - list with names as one char strings and 1 or 2 as turn order
      # @param: grid_size - numeric vector with two equivalent, odd integers
      # Returns: None
      
      self$EMPTY_CHAR <- empty_char
      self$CHARS <- chars
      self$GRID_SIZE <- grid_size
      
      private$setup()
    },

    start = function(){
      
      # Main game loop start.
      # @param: None
      # Returns: None
      
      while (is.null(private$winner)) {
        if (private$current_player == private$player_2 & isTRUE(private$computer)) {
          private$comp_move()
        } else {
          private$place_symbol(private$current_player)
        }
        
        # Set next player
        private$current_player <- names(self$CHARS)[private$current_player != names(self$CHARS)]
      
        private$.turn <- private$.turn + 1
        
        private$win_status <- private$check_win(self$grid)
        private$tie_status <- private$check_tie()
        
        if (isTRUE(private$win_status$Status)) {
          private$winner <- private$win_status$Chr
          
          # Show final board.
          print("\n")
          private$print_grid()
          
          if (isTRUE(private$computer & private$winner != private$player_1)) {
            cat("You lose...\n")
          } else {
            cat(sprintf("Congrats (%s)!\n", private$winner))
          }
        } else if (isTRUE(private$tie_status)) {
          private$winner <- "Tie"
          cat("Tie game!\n")
        }
      
        if (isTRUE(private$win_status$Status) | isTRUE(private$tie_status)) {
          cat(sprintf("Game lasted %s turns.\n", private$turns))
          cat("New game? Enter 'c' to continue.\n")
          new_game <- readLines("stdin", 1)
          if (new_game == "c") {
            private$setup()
            self$start()
          }
        }
      }
    }
  ),
  
  active = list(
    EMPTY_CHAR = function(value) {
      if (missing(value)) {
        private$.EMPTY_CHAR
      } else {
        if (!is.null(private$.EMPTY_CHAR)) {
          stop("Empty character already set.")
        } else {
          stopifnot(is.character(value))
          private$.EMPTY_CHAR <- value
        }
      }
    },
    CHARS = function(value) {
      if (missing(value)) {
        private$.CHARS
      } else {
        if (!is.null(private$.CHARS)) {
          stop("Characters already set.")
        } else {
          # Stop if not a list of one char vals of length 2
          stopifnot(is.list(value), 
                    length(value) == 2, 
                    all(sapply(names(value), length) == 1))
          # Stop if values not order nums 1 and 2. Checks if duplicates as well.
          if (!all(c(1, 2) %in% as.vector(unlist(value)))) {
            stop("Invalid order for characters. 1 or 2 only.")
          }
          private$.CHARS <- value
        }
      }
    },
    
    # Reverse characters list into named chr vector. 
    # Used to print grid.
    INV_CHARS = function(value) {
      if (missing(value)) {
        inv_chars <- names(self$CHARS)
        names(inv_chars) <- unlist(self$CHARS)
        inv_chars
      } else {
        stop("INV_CHARS is read-only.")
      }
    }
    ,
    GRID_SIZE = function(value) {
      if (missing(value)) {
        private$.GRID_SIZE
      } else {
        if (!is.null(private$.GRID_SIZE)) {
          stop("Grid size already set.")
        } else {
          # Stop if not a numeric vector of length 2
          stopifnot(is.numeric(value), length(value) == 2)
          # Stop if values not identical and not odd.
          if (!all(value[1] == value & value %% 2 != 0)) {
            stop("Invalid grid size. Must identical, odd values.")
          }
          private$.GRID_SIZE <- value
        }
      }
    },
    
    grid = function(value) {
      if (missing(value)) {
        private$.grid
      } else {
        stop("Grid is read-only.")
      }
    }, 
    
    turn = function(value) {
      if (missing(value)) {
        private$.turn
      } else {
        stop("Turn is read-only.")
      }
    }
  ),
  
  private = list (
    # Active but read-only. Constants set at init and cannot be changed after.
    .EMPTY_CHAR = NULL,
    .CHARS = NULL,
    .GRID_SIZE = NULL,
    .grid = NULL,
    .turn = 0,
    
    # Private fields
    computer = FALSE,
    player_1 = NULL,
    player_2 = NULL,
    current_player = NULL,
    
    winner = NULL,
    win_status = FALSE,
    tie_status = FALSE,
    
    setup = function(){
      
      # Setup game settings.
      # @param: None
      # Returns: None
      
      cat("Welcome to Tic-Tac-Toe!\n")
      cat("Play against a computer? (y/n)\n")
      comp_prompt <- readLines("stdin", 1)
      if (comp_prompt == "y") {
        private$computer <- TRUE
      } else {
        private$computer <- FALSE
      }
      
      private$player_1 <- private$intro()
      private$player_2 <- names(self$CHARS)[names(self$CHARS) != private$player_1]
      
      # Set order of turns.
      private$current_player <- ifelse(self$CHARS[[private$player_1]] == 1, private$player_1, private$player_2)
      
      # Reset winner if reset game.
      private$winner <- NULL

      # Create grid of NAs
      private$.grid <- matrix(data=NA, nrow=self$GRID_SIZE[1], ncol=self$GRID_SIZE[2])
    },
    
    intro = function(){
      
      # Game introduction prompt that sets player character.
      # @param: None
      # Returns: symbol [str]
      
      cat(paste0("Play as ", paste(names(self$CHARS), collapse = " or "), "?\n"))
      symbol <- readLines("stdin", 1)
      if (!symbol %in% names(self$CHARS)) {
        cat("Invalid character. Try again.\n")
        return(private$intro())
      } else {
        return(symbol)
      }
    },
    
    print_grid = function() {
      
      # Converts matrix grid to dataframe and prints.
      # @param: None
      # Returns: None
      
      chr_grid <- self$grid
      # Replace grid NAs with empty_char
      chr_grid[is.na(chr_grid)] <- self$EMPTY_CHAR
      
      # Replace 1s and 2s with given characters
      inv_chars <- self$INV_CHARS
      for (n in names(inv_chars)) {
        chr_grid[chr_grid == n] <- inv_chars[n]
      }
      
      # Convert to dataframe for prettier printing
      df_grid <- as.data.frame(chr_grid)
      names(df_grid) <- seq(self$GRID_SIZE[2])
      
      print(df_grid)
    },
    
    place_symbol = function(char) {
      
      # Player prompt to place symbol on grid. Also validates input.
      # @param: char - Player character [str]
      # Returns: None
      
      cat(sprintf("\nPlayer (%s):\n", private$current_player))
      private$print_grid()
      cat("Select a row:\n")
      row <- readLines("stdin", 1)
      cat("Select a col:\n")
      col <- readLines("stdin", 1)
      
      max_rows <- dim(self$grid)[1]
      max_cols <- dim(self$grid)[2]
      
      # convert to int. as.int coerces non-ints to NA
      row <- as.integer(row)
      col <- as.integer(col)
      
      if (!is.integer(row) | !is.integer(col) | is.na(row) | is.na(col)) { 
        cat("Invalid position (Not a valid number). Try again.")
        private$place_symbol(char)
        return()
      }
      
      if (col <= max_cols & row <= max_rows) {
        if (is.na(self$grid[row, col])) {
          private$.grid[row, col] <- self$CHARS[[char]]
        } else {
          cat("Invalid position (Already filled). Try again.")
          private$place_symbol(char)
          return()
        }
      } else {
        cat("Invalid position (Outside of bounds). Try again.")
        private$place_symbol(char)
        return()
      }
      
    },
    
    comp_move = function() {
      open_pos <- which(is.na(self$grid), arr.ind = TRUE)
      open_pos <- as.data.frame(open_pos)
      
      max_row <- dim(self$grid)[1]
      max_col <- dim(self$grid)[2]
      sel_row <- NULL
      sel_col <- NULL
      
      for (i in seq(dim(open_pos)[1])) {
        coord <- open_pos[i, ]
        row_ind <- coord[["row"]]
        col_ind <- coord[["col"]]
        
        grid_copy_p1 <- self$grid
        grid_copy_p2 <- self$grid
        
        # Create grid one move ahead for each player.
        grid_copy_p1[row_ind, col_ind] <- self$CHARS[[private$player_1]]
        grid_copy_p2[row_ind, col_ind] <- self$CHARS[[private$player_2]]
        
        # Once winning move found, immediately break from loop to prevent overwrite.
        # First, check available moves and chose move that would win.
        if (isTRUE(private$check_win(grid_copy_p1)$Status)) {
          sel_row <- row_ind
          sel_col <- col_ind
          break
          
        # Otherwise, check available moves and block any opponent winning moves.
        } else if (isTRUE(private$check_win(grid_copy_p2)$Status)) {
          sel_row <- row_ind
          sel_col <- col_ind
        }
      }
      
      # If no winning moves from either player, choose center, corner, or random.
      if (is.null(sel_row) & is.null(sel_col)) {
        open_corners <- open_pos[open_pos$row %in% c(1, max_row) & open_pos$col %in% c(1,max_col), ]
        
        # hacky soln to find center
        center <- c(median(seq(max_row)), median(seq(max_col)))
        
        # If center is open, else if corner free, then random.
        if (any(open_pos["row"] == center[1] & open_pos["col"] == center[2])) {
          sel_row <- center[1]
          sel_col <- center[2]
        } else if (length(open_corners) > 0) {
          # randomly select a corner
          corner_coord <- open_corners[sample(dim(open_corners)[1], 1), ]
          sel_row <- corner_coord[["row"]]
          sel_col <- corner_coord[["col"]]
        } else {
          random_pos <- open_pos[sample(dim(open_pos)[1], 1), ]
          sel_row <- random_pos[["row"]]
          sel_col <- random_pos[["col"]]
        }
        
      }

      # Place on grid.
      private$.grid[sel_row, sel_col] <- self$CHARS[[private$current_player]]
      
    },
    
    check_tie = function() {
      
      # Check if a tie based on if any NAs left on grid.
      # @param: None
      # Returns: True or False [bool]
      
      if (!any(is.na(self$grid))) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    },
    
    # Couldn't implement python version.
    # Check if uniq values in row/col is 1 (ie. all the same)
    # https://stackoverflow.com/a/35567910
    check_win = function(grid){
      
      # Checks win through a row or column.
      # @param: grid - matrix of grid
      # Returns: status - list with Status [bool] and player character or "None" [str]
      
      status <- list("Status" = FALSE, "Chr" = NA)
      
      win_rows <- grid[apply(grid, 1, function(x) length(unique(x)) == 1 & !any(is.na(x))), ]
      win_cols <- grid[, apply(grid, 2, function(x) length(unique(x)) == 1 & !any(is.na(x)))]
      
      diag_1 <- diag(grid)
      # flip grid by cols
      diag_2 <- diag(grid[,ncol(grid):1])

      if (length(win_rows) != 0) {
        status["Chr"] <- self$INV_CHARS[as.character(win_rows[1])]
      } else if(length(win_cols) != 0) {
        status["Chr"] <- self$INV_CHARS[as.character(win_cols[1])]
      } 
      
      # all() will return NAs if any NAs in vector. Remove and second check with correct diagonal length.
      if (all(diag_1 == diag_1[1], na.rm = TRUE) & length(diag_1[!is.na(diag_1)]) == self$GRID_SIZE[1]) {
        status["Chr"] <- self$INV_CHARS[as.character(diag_1[1])]
      } else if (all(diag_2 == diag_2[1], na.rm = TRUE) & length(diag_2[!is.na(diag_2)]) == self$GRID_SIZE[1]) {
        status["Chr"] <- self$INV_CHARS[as.character(diag_2[1])]
      }
      
      # Set win status to TRUE if character set.
      if (!is.na(status[["Chr"]])) {
        status["Status"] <- TRUE
      }
      
      return(status)
    }
    
  )

)

game1 <- TicTacToe$new("*", list("X" = 1, "O" = 2), c(3, 3))
game1$start()
