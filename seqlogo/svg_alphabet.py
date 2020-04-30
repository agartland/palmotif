"""Most are just paths with fill.
Second and third paths are white e.g. inside of A.
O and Q require circles and S requires stroke only.

Credit for the SVG to creators of logojs:
https://github.com/weng-lab/logojs-package"""

__all__ = ['svg_letter',
           'add_letter']

alphabet_paths = {'A':['M 0 100 L 33 0 L 66 0 L 100 100 L 75 100 L 66 75 L 33 75 L 25 100 L 0 100',
                        'M 41 55 L 50 25 L 58 55 L 41 55'],
                  'B':['M 0 0 L 80 0 C 105 0 105 50 80 50 C 105 50 105 100 80 100 L 00 100 L 0 0',
                       'M 20 15 L 70 15 C 80 15 80 35 70 35 L 20 35 L 20 15',
                       'M 20 65 L 70 65 C 80 65 80 85 70 85 L 20 85 L 20 65'],
                  'C':['M 100 28 C 100 -13 0 -13 0 50 C 0 113 100 113 100 72 L 75 72 C 75 90 30 90 30 50 C 30 10 75 10 75 28 L 100 28'],
                  'D':['M 0 0 L 60 0 C 110 0 110 100 60 100 L 0 100 L 0 0',
                       'M 20 15 L 40 15 C 85 15 85 85 40 85 L 20 85 L 20 15'],
                  'E':['M 0 0 L 100 0 L 100 20 L 20 20 L 20 40 L 90 40 L 90 60 L 20 60 L 20 80 L 100 80 L 100 100 L 0 100 L 0 0'],
                  'F':['M 0 0 L 100 0 L 100 20 L 20 20 L 20 40 L 80 40 L 80 60 L 20 60 L 20 100 L 0 100 L 0 0'],
                  'G':['M 100 28 C 100 -13 0 -13 0 50 C 0 113 100 113 100 72 L 100 48 L 55 48 L 55 72 L 75 72 C 75 90 30 90 30 50 C 30 10 75 5 75 28 L 100 28'],
                  'H':['M 0 0 L 20 0 L 20 40 L 80 40 L 80 0 L 100 0 L 100 100 L 80 100 L 80 60 L 20 60 L 20 100 L 0 100 L 0 0'],
                  'I':['M 40 0 L 60 0 L 60 100 L 40 100 L 40 0'],
                  'J':['M 0 60 C 0 111 100 111 100 60 L 100 0 L 75 0 L 75 60 C 80 90 20 90 25 60'],
                  'K':['M 0 0 L 20 0 L 20 40 L 75 0 L 100 0 L 50 50 L 100 100 L 75 100 L 30 65 L 20 75 L 20 100 L 0 100 L 0 0'],
                  'L':['M 0 0 L 0 100 L 100 100 L 100 80 L 20 80 L 20 0 L 0 0'],
                  'M':['M 0 0 L 20 0 L 50 35 L 80 0 L 100 0 L 100 100 L 80 100 L 80 30 L 50 65 L 20 30 L 20 100 L 0 100 L 0 0'],
                  'N':['M 0 100 L 0 0 L 20 0 L 80 75 L 80 0 L 100 0 L 100 100 L 80 100 L 20 25 L 20 100 L 0 100'],
                  'O':['cx="50" cy="50" r="50"', 'cx="50" cy="50" r="32"'], # requires circles
                  'P':['M 0 0 L 80 0 C 105 0 105 50 80 50 L 20 50 L 20 100 L 0 100 L 0 0',
                       'M 20 15 L 70 15 C 80 15 80 35 70 35 L 20 35 L 20 15'],
                  'Q':['M 85 100 L 55 70 L 70 55 L 100 85 L 85 100', 'cx="50" cy="50" r="50"', 'cx="50" cy="50" r="32"'], # O with slash
                  'R':['M 0 0 L 80 0 C 105 0 105 50 80 50 C 100 50 100 70 100 70 L 100 100 L 80 100 L 80 80 C 80 80 80 60 50 60 L 20 60 L 20 100 L 0 100 L 0 0',
                        'M 20 15 L 70 15 C 80 15 80 35 70 35 L 20 35 L 20 15'],
                  'S':['M92 26 A43 20 0 1 0 43 46 A42 23 0 1 1 9 68'], #stroke only
                  'T':['M 0 0 L 0 20 L 35 20 L 35 100 L 65 100 L 65 20 L 100 20 L 100 0 L 0 0'],
                  'U':['M 0 0 L 0 60 C 0 111 100 111 100 60 L 100 0 L 75 0 L 75 60 C 80 90 20 90 25 60 L 25 0 L 0 0'],
                  'V':['M 0 0 L 20 0 L 50 80 L 80 0 L 100 0 L 60 100 L 40 100 L 0 0'],
                  'W':['M 0 0 L 20 0 L 30 70 L 50 30 L 70 70 L 80 0 L 100 0 L 90 100 L 70 100 L 50 65 L 30 100 L 10 100 L 0 0'],
                  'X':['M 0 0 L 20 0 L 50 40 L 80 0 L 100 0 L 70 50 L 100 100 L 80 100 L 50 60 L 20 100 L 0 100 L 30 50 L 0 0'],
                  'Y':['M 0 0 L 20 0 L 50 45 L 80 0 L 100 0 L 60 60 L 60 100 L 40 100 L 40 60 L 0 0'],
                  'Z':['M 0 0 L 100 0 L 100 20 L 35 80 L 100 80 L 100 100 L 0 100 L 0 80 L 65 20 L 0 20 L 0 0'],
                  '-':['M 12.5 40 H 87.5 V 60 H 12.5 V 40']}

def svg_letter(a, color='black', background='white'):
    paths_only = 'ABCDEFGHIJKLMNPRTUVWXYZ-'
    if a in paths_only:
        out = [{'element':'path', 'd':alphabet_paths[a][0], 'fill':color}]
        for ap in alphabet_paths[a][1:]:
            out.append({'element':'path', 'd':ap, 'fill':background})
    elif a == 'O':
        out = [{'element':'circle', 'center':(50, 50), 'r':50, 'fill':color},
               {'element':'circle', 'center':(50, 50), 'r':32, 'fill':background}]
    elif a == 'Q':
        out = [{'element':'circle', 'center':(50, 50), 'r':50, 'fill':color},
               {'element':'circle', 'center':(50, 50), 'r':32, 'fill':background},
               {'element':'path', 'd':'M 85 100 L 55 70 L 70 55 L 100 85 L 85 100', 'fill':color}]
    elif a == 'S':
        out = [{'element':'path', 'd':'M92 26 A43 20 0 1 0 43 46 A42 23 0 1 1 9 68', 'stroke':color, 'fill':background, 'stroke-width':18}]
    return out

def add_letter(dwg, letter, group_id='', color='black', background='white', **attributes):
    if group_id is None:
        group_id = 'letter_%s' % letter
    elements = svg_letter(letter, color=color, background=background)
    group = dwg.add(dwg.g(id=group_id, **attributes))
    for e in elements:
        if e['element'] == 'path':
            tmp_attrib = {k:e[k] for k in e if not k == 'element'}
            group.add(dwg.path(**tmp_attrib))
        elif e['element'] == 'circle':
            tmp_attrib = {k:e[k] for k in e if not k == 'element'}
            group.add(dwg.circle(**tmp_attrib))
    return group


