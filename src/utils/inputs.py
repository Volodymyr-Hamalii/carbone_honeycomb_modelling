class Inputs:
    @staticmethod
    def bool_input(to_set: bool, default_value: bool, text: str) -> bool:
        if to_set is False:
            return default_value

        default_value_str: str = "+" if default_value else "-"
        message: str = f"{text} (print '+', '-' or nothing; by default is '{default_value_str}'): "
        result_str: str = input(message) or default_value_str
        return result_str == "+"

    @staticmethod
    def text_input(to_set: bool, default_value: str, text: str):
        if to_set is False:
            return default_value

        message: str = f"{text} (by default is '{default_value}'): "
        return input(message) or default_value
