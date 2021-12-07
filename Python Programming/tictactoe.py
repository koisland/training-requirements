import numpy as np
import pandas as pd


class TicTacToe:
    EMPTY_CHAR = "*"
    # 1 determines turn order as well
    CHARS = {"X": 1, "O": 0}
    # Both row and cols need to be odd.
    GRID_SIZE = (5, 5)
    INV_CHARS = {v: k for k, v in CHARS.items()}

    def __init__(self):
        self._setup()

    def _setup(self):
        print("Welcome to Tic-Tac-Toe!\n")
        if input("Play against a computer? (y/n)\n") == "y":
            self.computer = True
        else:
            self.computer = False

        self.player_1 = self.intro()
        self.player_2 = [char for char in self.CHARS.keys() if char != self.player_1][0]
        # init grid with nans
        self.grid = np.empty(self.GRID_SIZE) * np.nan
        self.turn = 0
        self.curr_player = (self.player_1 if self.CHARS[self.player_1] == 1 else self.player_2)
        self.winner = None

    def intro(self):
        symbol = input(f"Play as {' or '.join(self.CHARS.keys())}?\n")
        if symbol not in self.CHARS.keys():
            print("Invalid character. Try again.\n")
            return self.intro()
        else:
            return symbol

    @property
    def grid(self):
        return self._grid

    @grid.setter
    def grid(self, grid):
        self._grid = grid

    def check_win_diag(self, grid):
        # get diagonal chrs
        diag_1 = np.diag(grid)
        diag_2 = np.diag(np.flip(grid, axis=1))

        # if all elements in the diagonal are equal to the first element, return winner.
        if np.all(diag_1 == diag_1[1]):
            return self.INV_CHARS.get(diag_1[1])
        elif np.all(diag_2 == diag_2[1]):
            return self.INV_CHARS.get(diag_2[1])
        else:
            # otherwise, no winner and return None.
            return

    def check_win_row_col(self, grid):
        # get first values in row/col and compare to rest of matrix.
        # if all equal in row/col, then return logical matrix with winning row/col
        cols = np.all(grid == grid[0, :], axis=0)
        # transpose grid and compare to first value of row
        rows = np.all(grid.T == grid[:, 0], axis=0)

        # nonzero gives indices where not zero ie. True (1)
        # https://stackoverflow.com/questions/16094563/numpy-get-index-where-value-is-true
        win_row = np.nonzero(rows)[0]
        win_col = np.nonzero(cols)[0]

        # slice grid to get winning row/col. if not empty, return first value (winner).
        # make sure to use chrs not num as may eval as false (0)
        if grid[:, win_col].size > 0:
            return self.INV_CHARS.get(grid[:, win_col][0, 0])
        elif grid[win_row, :].size > 0:
            return self.INV_CHARS.get(grid[win_row, :][0, 0])
        else:
            return

    def check_tie(self):
        if not np.any(np.isnan(self.grid)):
            return True

    def comp_move(self):
        # open positions as a list of tuples w/row and col
        open_pos = list(zip(*np.nonzero(np.isnan(self.grid))))
        # max row and col indices for corners.
        max_row_ind, max_col_ind = np.array(self.grid.shape) - 1
        # selected row and col
        sel_row, sel_col = None, None

        for row, col in open_pos:
            grid_copy_p1 = self.grid.copy()
            grid_copy_p2 = self.grid.copy()
            grid_copy_p1[row][col] = self.CHARS[self.player_1]
            grid_copy_p2[row][col] = self.CHARS[self.player_2]

            # Check available moves and chose move that would win.
            if self.check_win_row_col(grid_copy_p2) or self.check_win_diag(grid_copy_p2):
                sel_row, sel_col = row, col
            # Check available moves and block any opponent winning moves.
            elif self.check_win_row_col(grid_copy_p1) or self.check_win_diag(grid_copy_p1):
                sel_row, sel_col = row, col

        # if no winning moves from either player, choose center, corner, or random in that order.
        if sel_row is None and sel_col is None:
            open_corners = [(r, c) for r, c in open_pos if r in (0, max_row_ind) and c in (0, max_col_ind)]

            # Check if center is open.
            if (max_row_ind / 2, max_col_ind / 2) in open_pos:
                sel_row, sel_col = (int(max_row_ind / 2), int(max_col_ind / 2))

            # Check if a corner is open.
            elif len(open_corners) > 0:
                corner_indice = np.random.choice(range(len(open_corners)))
                sel_row, sel_col = open_corners[corner_indice]

            # Otherwise, choose random pos.
            else:
                pos_indice = np.random.choice(range(len(open_pos)))
                sel_row, sel_col = open_pos[pos_indice]

        # make grid placement
        self.grid[sel_row][sel_col] = self.CHARS[self.curr_player]

    def place_symbol(self, char):
        print(f"\nPlayer ({self.curr_player}):")
        self.print_grid()
        row = input("Select a row:\n")
        col = input("Select a column:\n")

        max_row_ind, max_col_ind = np.array(self.grid.shape) - 1

        # check if str is number. isnumeric and isdigit only ork on unicode
        try:
            row = int(row)
            col = int(col)
        except ValueError:
            print("Invalid position (Not a valid number). Try again.")
            self.place_symbol(char)

        # check if below or equal to max ind
        if int(col) <= max_col_ind and int(row) <= max_row_ind:
            row_copy = self.grid[row]
            # if position is empty.
            if np.isnan(row_copy[col]):
                row_copy[col] = self.CHARS[char]
            else:
                print("Invalid position (Already filled). Try again.")
                self.place_symbol(char)

            self.grid[row] = row_copy
        else:
            print("Invalid position (Outside of bounds). Try again.")
            self.place_symbol(char)

    def print_grid(self):
        # cvt numpy matrix to pd dataframe to print
        grid_df = pd.DataFrame(self.grid)
        # fill NaNs with empty char
        grid_df.fillna(self.EMPTY_CHAR, inplace=True)

        # replace 0 or 1 in df with chars.
        for k, v in self.INV_CHARS.items():
            grid_df.replace(k, v, inplace=True)

        print(grid_df, "\n")

    def start(self):
        while self.winner is None:
            if self.curr_player == self.player_2 and self.computer == True:
                self.comp_move()
            else:
                self.place_symbol(self.curr_player)

            # bad. find a better way to figure out
            self.curr_player = [char for char in self.CHARS.keys() if char != self.curr_player][0]
            self.turn += 1

            win_status = self.check_win_row_col(self.grid) or self.check_win_diag(self.grid)
            tie_status = self.check_tie()

            # check winner in row, col, or diags or tie.
            if win_status:
                self.winner = win_status
                # if playing against a computer and user is not the winner.
                if self.computer and self.winner != self.player_1:
                    print("You lose...")
                else:
                    print(f"Congrats ({self.winner})!")
            elif tie_status:
                self.winner = "Tie"
                print("Tie game!")

            # reset game.
            if win_status or tie_status:
                print(f"Game lasted {self.turn} turns.\n")
                new_game = input("New game? Enter 'c' to continue.\n")
                if new_game == "c":
                    self._setup()


def main():
    ttt_game = TicTacToe()
    ttt_game.start()


if __name__ == "__main__":
    main()
