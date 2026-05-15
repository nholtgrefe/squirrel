
## Global variable indicating whether validation (i.e. checks whether input is correct
## form etc.) is on for the method. Saves time if False.

VALIDATE = False

def set_validate(value: bool):
    global VALIDATE
    VALIDATE = value

def validate() -> bool:
    return VALIDATE