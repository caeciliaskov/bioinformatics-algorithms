import time

S = "hpphhhphhh"

def evens_odds(S):
   evens_odds = []
   pointer = 0
   for i in S:
      if i == 'h' and pointer % 2 == 0:
         evens_odds.append('e')
      if i == 'h' and pointer % 2 == 1:
         evens_odds.append('o')
      elif i != "h":
         evens_odds.append(" ")
      pointer += 1
   return evens_odds

print(evens_odds(S))

def match_evens_from_left(S):
   e_o_list = evens_odds(S)
   matches = []
   i = 0
   j = len(e_o_list)-1
   while i < 0.5 * len(e_o_list) and j > 0.5 * len(e_o_list):
      if e_o_list[i] == 'e' and e_o_list[j] == 'o':
         matches.append((i, j))
         i += 1
         j -= 1
      if e_o_list[i] == 'e' and e_o_list[j] != 'o':
         j -= 1
      if e_o_list[i] != 'e' and e_o_list[j] == 'o':
         i += 1
      if e_o_list[i] != 'e' and e_o_list[j] != 'o':
         i += 1
         j -= 1
   return matches
 
print(match_evens_from_left(S))

def match_odds_from_left(S):
   e_o_list = evens_odds(S)
   matches = []
   i = 0
   j = len(e_o_list)-1
   while i < 0.5 * len(e_o_list) and j > 0.5 * len(e_o_list):
      if e_o_list[i] == 'o' and e_o_list[j] == 'e':
         matches.append((i, j))
         i += 1
         j -= 1
      if e_o_list[i] == 'o' and e_o_list[j] != 'e':
         j -= 1
      if e_o_list[i] != 'o' and e_o_list[j] == 'e':
         i += 1
      if e_o_list[i] != 'o' and e_o_list[j] != 'e':
         i += 1
         j -= 1
   return matches

print(match_odds_from_left(S))

def max_match(S):
   if len(match_evens_from_left(S)) > len(match_odds_from_left(S)):
      max = match_evens_from_left(S)
   else: 
      max = match_odds_from_left(S)
   return max

print(max_match(S))

def hpfold(S):
   forward = max_match(S)
   fold = ""
   i = 0
   fold = fold + (forward[0][0]) * "f"
   while i in range(len(forward)-1):
      D = forward[i+1][0] - forward[i][0] 
      K = D - 4
      if D == 2:
         fold = fold + "ff"
      elif D > 2 and D % 2 == 0:
         fold = fold + "l" + (int(K/2) * "f") + "rr" + (int(K/2) * "f") + "l"
      i += 1
   D = forward[-1][1] - forward[-1][0]
   K = D // 2
   fold = fold + (int(K)*"f") + "rr" + (int(K-1)*"f")
   j = len(forward)-1
   while j in range(len(forward)-1, 0, -1):
      D = forward[j-1][1] - forward[j][1]
      K = D - 4
      if D == 2:
         fold = fold + "ff"
      elif D > 2 and D % 2 == 0:
         fold = fold + "l" + (int(K/2) * "f") + "rr" + (int(K/2) * "f") + "l"
      j -= 1
   fold = fold + (len(S)-1 - forward[0][1]) * "f"
   return fold

start = time.time()

print(hpfold(S))

end = time.time()

#print(end - start)

#print(max_match(S)[1][0])