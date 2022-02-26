# library(Matrix)
library(R6)


TicTacToe <- R6Class("TicTacToe",
  public = list(
    initialize = function(empty_char, chars, grid_size) {
      self$EMPTY_CHAR <- empty_char
      self$CHARS <- chars
      self$GRID_SIZE <- grid_size
      
      private$setup()
    },

    start = function(){
      while (is.null(private$winner)) {
        if (private$current_player == private$current_player & isTRUE(private$computer)) {
          private$comp_move()
        } else {
          private$place_symbol(private$current_player)
        }
        
        # Set next player
        private$current_player <- names(self$CHARS)[private$current_player != names(self$CHARS)]
      
        private$turn <- private$turn + 1
        
        private$win_status <- private$check_win_row_col() | private$check_win_diag() 
        private$tie_status <- private$check_tie()
        
        if (isTRUE(private$win_status)) {
          private$winner <- private$win_status
          
          if (isTRUE(private$computer & private$winner != private$player_1)) {
            cat("You lose...\n")
          } else {
            cat(sprintf("Congrats (%s)!\n", private$winner))
          }
        } else if (isTRUE(private$tie_status)) {
          private$winner <- "Tie"
          cat("Tie game!\n")
        }
        
        if (isTRUE(private$win_status) | isTRUE(private$tie_status)) {
          cat(sprintf("Game lastes %s turns.\n", private$turns))
          new_game <- readline("New game? Enter 'c' to continue.\n")
          if (new_game == "c") {
            private$setup()
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
          # Stop if not a list of length 2
          stopifnot(is.list(value), length(value) == 2)
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
      cat("Welcome to Tic-Tac-Toe!\n")
      comp_prompt <- readline("Play against a computer? (y/n)\n")
      if (comp_prompt == "y") {
        private$computer <- TRUE
      } else {
        private$computer <- FALSE
      }
      
      private$player_1 <- private$intro()
      private$player_2 <- names(self$CHARS)[private$player_1 != names(self$CHARS)]
      private$current_player <- ifelse(self$CHARS[private$player_1] == 2, private$player_1, private$player_2)
      
      # Create grid of NAs
      private$.grid <- matrix(data=NA, nrow=self$GRID_SIZE[1], ncol=self$GRID_SIZE[2])
    },
    
    intro = function(){
      symbol <- readline(paste0("Play as ", paste(names(self$CHARS), collapse = " or "), "?\n"))
      if (!symbol %in% names(self$CHARS)) {
        cat("Invalid character. Try again.\n")
        return(private$intro())
      } else {
        return(symbol)
      }
    },
    
    print_grid = function() {
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
    
    place_symbol = function() {
      
    },
    
    comp_move = function() {
      
    },
    
    check_tie = function() {
      
    },
    
    check_win_row_col = function(){
      
    },
    
    check_win_diag = function(){
      diag_1 <- diag(self$grid)
      # flip grid by cols
      diag_2 <- diag(self$grid[,ncol(self$grid):1])
      
      if (all(diag_1 == diag_1[1])) {
        
      } else if (all(diag_2 == diag[1])) {
        
      } else {
        return(FALSE)
      }
    }
    
  )

)

game1 <- TicTacToe$new("*", list("X" = 2, "O" = 1), c(3, 3))
game1$start()