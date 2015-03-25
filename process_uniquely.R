#reads lists from csv file and assemble into membership table for easier loading alter

# asssumed format is rows of csv file, each row is list name followed by members.
tmp = read.table("data_intermediate//uniqOverlaps_H3K4_2FE_JB03182015.csv", sep = "\t", stringsAsFactors = F)

all_members = character()
all_lists = list()
for (r in 1:nrow(tmp)) {
  str = tmp[r, ]
  str = strsplit(str, ",")[[1]]
  name = str[1]
  members = str[2:length(str)]
  keep = members != ""
  members = members[keep]
  all_members = union(all_members, members)
  all_lists[[length(all_lists) + 1]] = members
  names(all_lists)[length(all_lists)] = name
}

membershipTable = matrix(data = F, nrow = length(all_members), ncol = length(all_lists))
rownames(membershipTable) = all_members
colnames(membershipTable) = names(all_lists)
for (i in 1:length(all_lists)) {
  m = all_lists[[i]]
  membershipTable[m, i] = T
}
uniquely_K4_membership = membershipTable
save(uniquely_K4_membership, file = 'data_intermediate///uniquely_membership.save')
