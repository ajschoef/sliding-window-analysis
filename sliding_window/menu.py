

class Menu:

    def __init__(self):
        self.current_parameters
        self.menu_text = """================================\n
           MENU\n
               ================================\n
               1 - Set collinear neighbors filter (default = 0.1)\n
               2 - Set the upper limit of regularization parameter (default = 1e-3)\n
               3 - Set the lower limit of regularization parameter (default = 1e-7)\n
               4 - Set the max iterations for the iterative solver (default = 1000)\n
               5 - Set the number of cross-validation grid values (default = 20)
               6 - Fit model\n
               ================================\n
           Enter a choice and press enter: """

    def menu_options(self, option):
        if option == 1:
            input("Set collinear neighbors filter: ")
        elif option == 2:
            input("Set the upper limit of regularization parameter (default = 1e-3)")
        elif option == 3:
            input("Set the lower limit of regularization parameter (default = 1e-7)")
        elif option == 4:
            input("Set the max iterations for the iterative solver (default = 1000)")
        elif option == 5:
            input("Set the number of cross-validation grid values (default = 20)")

    def run_menu(self):
        option = None
        while True:
            option = input(self.menu_text)
            if int(option) == 6:
                print("Fitting model...")
                break
