import  torch
from    torch import nn
from    torch.nn import functional as F
from    torch import optim
import  torch.utils.data as Data
from    sklearn.model_selection import KFold
import  numpy as np
import  pickle as pkl
from model import Transformer
from sklearn.metrics import classification_report
from sklearn.metrics import precision_score, recall_score

transformer_fn = 'transformer_input/'
pcs2idx = pkl.load(open(transformer_fn+'pc2wordsid.dict', 'rb'))
num_pcs = len(set(pcs2idx.keys()))


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
torch.cuda.set_device(0)
src_pad_idx = 0
src_vocab_size = num_pcs+1



def reset_model():
    model = Transformer(
                src_vocab_size, 
                src_pad_idx, 
                device=device, 
                max_length=300, 
                dropout=0.1
    ).to(device)
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    loss_func = nn.BCEWithLogitsLoss()
    return model, optimizer, loss_func



def return_batch(train_sentence, label, flag):
    X_train = torch.from_numpy(train_sentence).to(device)
    y_train = torch.from_numpy(label).float().to(device)
    train_dataset = Data.TensorDataset(X_train, y_train)
    training_loader = Data.DataLoader(
        dataset=train_dataset,    
        batch_size=200,
        shuffle=flag,               
        num_workers=0,              
    )
    return training_loader



def return_tensor(var, device):
    return torch.from_numpy(var).to(device)





def reject_prophage(all_pred, weight):
    all_pred = np.array(all_pred)
    all_pred[weight < 0.3] = 0
    return all_pred



# training with short contigs 
model, optimizer, loss_func = reset_model()
try:
    pretrained_dict=torch.load('database/transformer.pth', map_location=device)
    model.load_state_dict(pretrained_dict)
except:
    print('cannot find pre-trained model')
    exit()

####################################################################################
##########################    train with contigs    ################################
####################################################################################


sentence   = open(transformer_fn + 'sentence.feat', 'wb')
id2contig  = open(transformer_fn + 'sentence_id2contig.dict', 'wb')
proportion = open(transformer_fn + 'sentence_proportion.feat', 'wb')


all_pred = []
all_score = []
with torch.no_grad():
    for idx in range(0, len(sentence), 500):
        try:
            batch_x = sentence[idx: idx+500]
            weight  = proportion[idx: idx+500]
        except:
            batch_x = sentence[idx:]
            weight  = proportion[idx:]
        batch_x = return_tensor(batch_x, device).int()
        logit = model(batch_x)
        logit = torch.sigmoid(logit.squeeze(1))
        logit = reject_prophage(logit, weight)
        pred  = ['phage' if item > 0.5 else 'non-phage' for item in logit]
        all_pred += pred
        all_score += [float('{:.3f}'.format(i)) for i in logit]


pred_csv = pd.DataFrame({"Contig":id2contig.values(), "Pred":all_pred, "Score":all_score})
pred_csv.to_csv('final_prediction.csv', index = False)




















