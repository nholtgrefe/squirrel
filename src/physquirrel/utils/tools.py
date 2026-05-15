import string, random

############################################################

def id_generator(nr_ids=1, string_len=10):
    """Generate a list of 'nr_ids' random id-strings of size 'string_len'. If
    'nr_ids = 1', a single string (not in a list) is returned."""
    chars = string.ascii_uppercase + string.ascii_lowercase + string.digits
    if nr_ids == 1:
        return ''.join(random.choice(chars) for _ in range(string_len))
    else:
        return [id_generator(string_len=string_len) for _ in range(nr_ids)]
