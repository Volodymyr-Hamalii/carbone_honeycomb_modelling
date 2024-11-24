import sys
from src import AppActions, CommandLineArgsHandler, Constants, Logger


logger = Logger("Main")


def run_action() -> None:
    action: str = CommandLineArgsHandler.get_str_arg(1, default=Constants.DEFAULT_ACTION)
    structure_folder: str = CommandLineArgsHandler.get_str_arg(2, default=Constants.DEFAULT_STRUCTURE_FOLDER)
    to_set: bool = CommandLineArgsHandler.get_bool_arg(3, to_compare="set", default=False)
    args: list[str] = sys.argv[4:]

    if action == "help":
        AppActions.help(structure_folder)
        return

    try:
        logger.info(f"Run '{action}' for '{structure_folder}' structure.")

        # Get the method from the AppActions class
        method = getattr(AppActions, action)

        if args:
            # Call the method with the structure_folder and other parameters
            method(structure_folder, to_set, *args)
        else:
            # Call the method with the structure_folder sparameter
            method(structure_folder, to_set)

    except AttributeError:
        logger.error(f"Action '{action}' is not available.")
        AppActions.help(structure_folder)

    except TypeError as e:
        logger.error(f"Parameters error:", e)
        AppActions.help(structure_folder)


if __name__ == "__main__":
    run_action()
