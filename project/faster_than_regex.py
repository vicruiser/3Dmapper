import re
with open('/usr/share/dict/american-english') as wordbook:
    banned_words = set(word.strip().lower() for word in wordbook)


def delete_banned_words(matchobj):
    word = matchobj.group(0)
    if word.lower() in banned_words:
        return ""
    else:
        return word

sentences = ["I'm eric. Welcome here!", "Another boring sentence.",
             "GiraffeElephantBoat", "sfgsdg sdwerha aswertwe"] * 250000

word_pattern = re.compile('\w+')

for sentence in sentences:
    sentence = word_pattern.sub(delete_banned_words, sentence)