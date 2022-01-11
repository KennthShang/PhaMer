import  torch
from    torch import nn
from    torch.nn import functional as F
from    torch import optim
import  torch.utils.data as Data
from    sklearn.model_selection import KFold
import  numpy as np
import  pickle as pkl
from model import Transformer
from return_attention import Transformer
from sklearn.metrics import classification_report
from sklearn.metrics import precision_score, recall_score


pcs2idx = pkl.load(open('output/pc2wordsid.dict', 'rb'))
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
    all_pred[weight < 0.2] = 0
    return all_pred

# training with short contigs 
model, optimizer, loss_func = reset_model()
try:
    pretrained_dict=torch.load('database/transformer.pth', map_location=device)
    model.load_state_dict(pretrained_dict)
except:
    print('cannot find pre-trained model')

####################################################################################
##########################    train with contigs    ################################
####################################################################################


train_phage           = pkl.load(open('dataset/2018.feat', 'rb'))
train_contig          = pkl.load(open('contigs/2018_10k.feat', 'rb'))
train_phage_weight    = pkl.load(open('dataset/2018_proportion.feat', 'rb'))
train_contig_weight   = pkl.load(open('contigs/2018_10k_proportion.feat', 'rb'))
train_bacteria        = pkl.load(open('dataset/train_bacteria.feat', 'rb'))
train_bacteria_weight = pkl.load(open('dataset/train_bacteria_proportion.feat', 'rb'))
test_phage            = pkl.load(open('contigs/2019_5k.feat', 'rb'))
test_phage_weight     = pkl.load(open('contigs/2019_5k_proportion.feat', 'rb'))
test_bacteria         = pkl.load(open('dataset/test_bacteria.feat', 'rb'))
test_bacteria_weight  = pkl.load(open('dataset/test_bacteria_proportion.feat', 'rb'))



train_label           = np.concatenate((np.ones(train_phage.shape[0]), np.ones(train_contig.shape[0]), np.zeros(train_bacteria.shape[0]), np.zeros(train_bacteria.shape[0])))
train_sentense        = np.vstack((train_phage, train_contig, train_bacteria, train_bacteria))
train_weight          = np.concatenate((train_phage_weight, train_contig_weight, train_bacteria_weight, train_bacteria_weight))
train                 = train_sentense

test_label            = np.concatenate((np.ones(test_phage.shape[0]), np.zeros(test_bacteria.shape[0])))
test_sentense         = np.vstack((test_phage, test_bacteria))
test_weight           = np.concatenate((test_phage_weight, test_bacteria_weight))
test                  = test_sentense



training_loader = return_batch(train, train_label, flag = True)
test_loader = return_batch(test, test_label, flag = False)
for epoch in range(10):
    _ = model.train()
    for step, (batch_x, batch_y) in enumerate(training_loader): 
        sentense = batch_x.int()
        prediction = model(sentense)
        loss = loss_func(prediction.squeeze(1), batch_y)
        optimizer.zero_grad()
        loss.backward() 
        optimizer.step()
    _ = model.eval()
    with torch.no_grad():
        all_pred = []
        for step, (batch_x, batch_y) in enumerate(test_loader): 
            sentense = batch_x.int()
            logit = model(sentense)
            logit  = torch.sigmoid(logit.squeeze(1))
            pred  = [1 if item > 0.5 else 0 for item in logit]
            all_pred += pred
        precision = precision_score(test_label, all_pred)
        recall = recall_score(test_label, all_pred)
    print(f'epoch no. {epoch} || precision: {precision} || recall: {recall}')

