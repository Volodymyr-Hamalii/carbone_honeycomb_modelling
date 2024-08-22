from .logger import Logger

logger = Logger(__name__)


class Inputs:
    @classmethod
    def bool_input(cls, to_set: bool, default_value: bool, text: str) -> bool:
        if to_set is False:
            return default_value

        available_values: list[str] = ["+", "-"]

        default_value_str: str = "+" if default_value else "-"
        message: str = f"{text} (print '+', '-' or nothing; by default is '{default_value_str}'): "
        result_str: str = input(message) or default_value_str

        if cls._validate_result(result_str, available_values):
            return result_str == "+"
        else:
            logger.warning("Wront input. Available values:", available_values)
            return cls.bool_input(to_set, default_value, text)

    @classmethod
    def text_input(cls, to_set: bool, default_value: str, text: str, available_values: list[str] = []):
        if to_set is False:
            return default_value

        message: str = f"{text} (by default is '{default_value}'): "
        result: str = input(message) or default_value

        if cls._validate_result(result, available_values):
            return result
        else:
            logger.warning("Wront input. Available values:", available_values)
            return cls.text_input(to_set, default_value, text, available_values)

    @staticmethod
    def _validate_result(result: str, available_values: list[str]) -> bool:
        if len(available_values) == 0 or result in available_values:
            return True
        return False
