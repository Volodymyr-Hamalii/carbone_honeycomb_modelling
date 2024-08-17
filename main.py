import sys
from src import Actions, Logger


logger = Logger(__name__)

DEFAULT_ACTION = "full_flow"
DEFAULT_STRUCTURE_FOLDER = "A1-7_h3"


def main() -> None:
    action: str = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] != "_" else DEFAULT_ACTION
    structure_folder: str = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] != "_" else DEFAULT_STRUCTURE_FOLDER
    to_set: bool = sys.argv[3] == "set" if len(sys.argv) > 3 and sys.argv[3] != "_" else False
    args: list[str] = sys.argv[4:]

    if action == "help":
        Actions.help(structure_folder)
        return

    try:
        logger.info(f"Run '{action}' for '{structure_folder}' structure.")

        # Get the method from the Actions class
        method = getattr(Actions, action)

        if args:
            # Call the method with the structure_folder and other parameters
            method(structure_folder, to_set, *args)
        else:
            # Call the method with the structure_folder sparameter
            method(structure_folder, to_set)

    except AttributeError:
        logger.error(f"Action '{action}' is not available.")
        Actions.help(structure_folder)

    except TypeError as e:
        logger.error(f"Parameters error:", e)
        Actions.help(structure_folder)


if __name__ == "__main__":
    main()
