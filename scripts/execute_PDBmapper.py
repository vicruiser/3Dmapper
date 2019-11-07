sub_df['range'] = sub_df.apply(
    lambda x: '-'.join(list(range(int(x['start']), int(x['end'])+1))), 1)
