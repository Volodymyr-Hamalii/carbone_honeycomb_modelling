import sys
from src import Actions, Logger


logger = Logger(__name__)

DEFAULT_ACTION = "full_flow"
DEFAULT_STRUCTURE_FOLDER = "A1-7_h3"


def main() -> None:
    action: str = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_ACTION
    structure_folder: str = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_STRUCTURE_FOLDER

    try:
        logger.info(f"Run '{action}' for '{structure_folder}' structure.")

        # Get the method from the Actions class
        method = getattr(Actions, action)

        # Call the method with the structure_folder parameter
        method(structure_folder)
    except AttributeError:
        logger.error(f"Action '{action}' is not available.")
        Actions.help(structure_folder)


if __name__ == "__main__":
    main()
