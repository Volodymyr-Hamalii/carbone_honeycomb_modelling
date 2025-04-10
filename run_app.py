from src import (
    Logger,
    AppGui,
)


logger = Logger("Main")


if __name__ == "__main__":
    app = AppGui()
    app.mainloop()
