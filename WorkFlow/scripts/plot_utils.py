import pandas as pd


def parse_indels(indel_str):
    positions_1bp = []
    positions = []
    if indel_str:
        entries = indel_str.split('|')
        for entry in entries:
            pos_type, _ = entry.split(':')
            pos, size = pos_type.split(
                'I') if 'I' in pos_type else pos_type.split('D')
            if size == "1":
                positions_1bp.append(int(pos))
            positions.append(int(pos))
    return positions, positions_1bp


def filter_reads(ins_str, del_str, lftS, lftE, rgtS, rgtE):
    """ 根据indel位置过滤reads """
    insertions, insertions_1bp = parse_indels(ins_str)
    deletions, deletions_1bp = parse_indels(del_str)
    # 合并插入和删除位置
    indel_positions = set(insertions + deletions)
    indel_positions_1bp = set(insertions_1bp + deletions_1bp)
    # 检查是否有indel在1-25或32-37范围内
    if any(pos in range(lftS-1,lftE) or pos in range(rgtS-1, rgtE) for pos in indel_positions_1bp):
        # 检查是否有indel在26-31范围内
        if any(pos in range(lftE, rgtS-1) for pos in indel_positions):
            return True  # 保留reads
        return False  # 删除reads
    return True


def check_wildtype(ins_str, del_str, mut_str):
    if ins_str == "" and del_str == "" and mut_str == "":
        return True
    return False


def parse_indels_size(indel_str):
    Size = []
    if indel_str:
        entries = indel_str.split('|')
        for entry in entries:
            pos_type, _ = entry.split(':')
            pos, size = pos_type.split(
                'I') if 'I' in pos_type else pos_type.split('D')
            Size.append(int(size))
    return sum(Size)


def parse_mutations_size(mut_str):
    entries = []
    if mut_str:
        entries = mut_str.split('|')
    return len(entries)


def make_frequency_table(fs,lftS, lftE, rgtS, rgtE, trapL):
    total_counts = 0
    Aligned_Sequence = []
    Reference_Sequence = []
    Unedited = []
    n_deleted = []
    n_inserted = []
    n_mutated = []
    Reads = []
    with open(fs, "r") as f:
        for line in f:
            if line.startswith("reference"):
                continue
            else:
                ref, query, count, ins, dels, mis = line.strip(
                    "\n").split("\t")
                Reference_Sequence.append(ref[:trapL])
                Aligned_Sequence.append(query[:trapL])
                n_inserted.append(parse_indels_size(ins))
                n_deleted.append(parse_indels_size(dels))
                n_mutated.append(parse_mutations_size(mis))
                total_counts += int(count)
                Une = True
                if filter_reads(ins, dels, lftS, lftE, rgtS, rgtE):
                    Une = False
                if check_wildtype(ins, dels, mis):
                    Une = True
                Unedited.append(Une)
                Reads.append(int(count))
    df = pd.DataFrame({"Aligned_Sequence": Aligned_Sequence, "Reference_Sequence": Reference_Sequence,
                      "Unedited": Unedited, "n_deleted": n_deleted, "n_inserted": n_inserted, "n_mutated": n_mutated, "Reads": Reads})
    df = df.assign(readP=df['Reads']/total_counts * 100)
    df.columns = ["Aligned_Sequence", "Reference_Sequence", "Unedited",
                  "n_deleted", "n_inserted", "n_mutated", "#Reads", "%Reads"]
    df = df.sort_values(by="#Reads", ascending=False, ignore_index=True)
    return df


def substitute_position(fs, refseq, lftS, lftE, rgtS, rgtE):
    # define a dictionary to store the refseq where the key is the position
    # and the value is the base
    refseq_dict = {}
    for i in range(len(refseq)):
        refseq_dict[i+1] = refseq[i]
    # define another dictionary to store the position and the base
    editing_dict = {}
    indel_dict = {}
    for i in range(len(refseq)):
        editing_dict[i+1] = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "-": 0}
    for i in range(len(refseq)):
        indel_dict[i+1] = {"Insertions": 0, "Insertions_Left": 0,
                           "Deletions": 0, "Substitutions": 0, "All_modifications": 0, "Total": 0}
    with open(fs, "r") as f:
        for line in f:
            if line.startswith("reference"):
                continue
            else:
                _, _, count, ins, dels, mis = line.strip("\n").split("\t")
                if filter_reads(ins, dels, lftS, lftE, rgtS, rgtE):
                    for _, v in indel_dict.items():
                        v['Total'] += int(count)
                    if mis != "":
                        mises = mis.split("|")
                        for mis in mises:
                            pos, mut = mis.split(":")
                            _, query_base = mut.split("2")
                            editing_dict[int(pos)][query_base] += int(count)
                            indel_dict[int(pos)]['Substitutions'] += int(count)
                    if ins != "":
                        ins_st = ins.split("|")
                        for ii in ins_st:
                            pos, mut = ii.split("I")
                            pos = int(pos)
                            if pos == 0:
                                pos = 1
                            indel_dict[pos]['Insertions_Left'] += int(count)
                            indel_dict[pos]['Insertions'] += int(count)
                    if dels != "":
                        del_st = dels.split("|")
                        for dd in del_st:
                            pos, mut = dd.split("D")
                            deln, _ = mut.split(":")
                            for i in range(int(deln)):
                                indel_dict[int(
                                    pos)+i]['Deletions'] += int(count)
                                editing_dict[int(pos)+i]['-'] += int(count)
    total_reads = indel_dict[1]['Total']
    for pp, v in editing_dict.items():
        v[refseq_dict[pp]] = total_reads - \
            (v['A'] + v['C'] + v['G'] + v['T'] + v['N'] + v['-'])
    for _, v in indel_dict.items():
        v['All_modifications'] = v['Insertions_Left'] + \
            v['Deletions'] + v['Substitutions']
    return editing_dict, indel_dict, total_reads
