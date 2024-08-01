import sys
from src import Actions, Logger


logger = Logger(__name__)

DEFAULT_ACTION = "full_flow"
DEFAULT_STRUCTURE_FOLDER = "A1-7_h3"


def main() -> None:
    if len(sys.argv) > 1:
        action: str = sys.argv[1]
        structure_folder: str = sys.argv[2]
    else:
        action = DEFAULT_ACTION
        structure_folder = DEFAULT_STRUCTURE_FOLDER

    try:
        # Get the method from the Actions class
        method = getattr(Actions, action)
        # Call the method with the structure_folder parameter
        method(structure_folder)
    except AttributeError:
        available_actions: list[str] = [
            func for func in dir(Actions)
            if callable(getattr(Actions, func))
            and not func.startswith("__")
        ]

        logger.error(f"Action '{action}' is not available. Available methods:", available_actions)


if __name__ == "__main__":
    main()
