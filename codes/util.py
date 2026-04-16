from torch.autograd import Variable
import torch.optim as optim
import torch
import torch.nn as nn
import numpy as np
import torch.nn.functional as f
import pandas as pd


if torch.cuda.is_available():
    device = "cuda"
elif torch.backends.mps.is_available() and torch.backends.mps.is_built():
    device = "mps"
else:
    device = "cpu"


class Encoder(nn.Module):
    def __init__(self, input_dim, hidden_dim, latent_dim):
        super(Encoder, self).__init__()
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_dim, latent_dim)

    def forward(self, x):
        x = self.relu(self.fc1(x))
        x = self.fc2(x)
        return x

class Decoder(nn.Module):
    def __init__(self, latent_dim, hidden_dim, output_dim):
        super(Decoder, self).__init__()
        self.fc1 = nn.Linear(latent_dim, hidden_dim)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_dim, output_dim)

    def forward(self, x):
        x = self.relu(self.fc1(x))
        x = self.fc2(x)
        return x

def RunAlgorithm(X, Y, anno, I, D=0, outpath = './', genelist = []):
    X_ori = X
    Y_ori = Y
    scName_full = X.columns.tolist()
    spName_full = Y.columns.tolist()
    Xallgenes = X.index.tolist()
    Yallgenes = Y.index.tolist()

    csX = X_ori.apply(lambda x: sum(x), axis=0)
    X_ori = X_ori / csX

    rsX = X.apply(lambda x: np.mean(x), axis=1)
    X = X.loc[rsX > 0,]
    rsX = X.apply(lambda x: np.mean(x), axis=1)
    X = X.loc[rsX < 10000 * np.median(rsX),]
    csX = X.apply(lambda x: sum(x), axis=0)
    X = X.loc[:, csX > np.median(csX) / 5]
    anno = anno[csX > np.median(csX) / 5]
    csX = X.apply(lambda x: sum(x), axis=0)
    X = X / csX

    rsY = Y.apply(lambda x: np.mean(x), axis=1)
    rsdY = Y.apply(lambda x: np.std(x), axis=1)
    covY = rsdY / (1 + rsY)
    Y = Y.loc[
        (rsY > np.median(rsY) / 5)
        & (rsY < np.mean(rsY) + 3 * np.std(rsY))
        & (covY > np.mean(covY) / 2),
    ]
    csY = Y.apply(lambda x: sum(x), axis=0)
    Y = Y / csY

    Xgenes = X._stat_axis.values.tolist()
    Ygenes = Y._stat_axis.values.tolist()
    commonGenes = list(set(list(set(Xgenes) & set(Ygenes))))
    if len(genelist) > 0:
        geneuse = list(set(commonGenes) & set(genelist))
    else:
        geneuse = commonGenes

    Xpart = X.loc[geneuse,]
    Ypart = Y.loc[geneuse,]

    scName_part = Xpart.columns.tolist()
    spName_part = Ypart.columns.tolist()

    Xdata = torch.tensor(Xpart.values).to(torch.float32).to(device)
    Ydata = torch.tensor(Ypart.values).to(torch.float32).to(device) 

    if D == 0:
        D_part = D.loc[spName_part, spName_part]
        I_part = I.loc[scName_part, scName_part]
        I_part_T = I_part.T
        I_part = (I_part + I_part_T + (I_part - I_part_T).abs()) / 2
        I_cut = np.percentile(I_part, 80)
        I_part[I_part < I_cut] = 0
        I_part[I_part < 0] = 0
        Idata = torch.tensor(I_part.values).to(torch.float32).to(device)
        Ddata = torch.tensor(D_part.values).to(torch.float32).to(device)
        Ddata = (Ddata - torch.mean(Ddata)) / torch.std(Ddata) + 2
        Ddata = torch.sigmoid(-Ddata / 2)

    input_dim = Xdata.shape[0]
    hidden_dim = 512
    latent_dim = 128
    encoder = Encoder(input_dim, hidden_dim, latent_dim).to(device)
    decoder = Decoder(latent_dim, hidden_dim, input_dim).to(device)

    num_epochs = 100
    learning_rate = 1e-3
    criterion = nn.MSELoss()
    params = list(encoder.parameters()) + list(decoder.parameters())
    optimizer = optim.Adam(params, lr=learning_rate)

    Xdata_train = Xdata.T 
    anno_arr = np.array(anno)
    alpha = 0.3

    for epoch in range(num_epochs):
        encoder.train()
        decoder.train()
        optimizer.zero_grad()
        latent = encoder(Xdata_train)
        reconstructed = decoder(latent)
        loss_recon = criterion(reconstructed, Xdata_train)

        latent_dist_loss = 0.0
        count = 0
        centers = []
        for label in np.unique(anno_arr):
            idx = np.where(anno_arr == label)[0]
            if len(idx) > 0:
                latent_group = latent[idx]
                center = latent_group.mean(dim=0, keepdim=True)
                centers.append(center)
                if len(idx) > 1:
                    dists = torch.norm(latent_group - center, dim=1) ** 2  
                    latent_dist_loss += dists.mean()
                    count += 1
        if len(centers) > 1:
            centers = torch.cat(centers, dim=0)
            center_dists = torch.cdist(centers, centers, p=2) ** 2 
            mask = torch.triu(torch.ones_like(center_dists), diagonal=1) > 0
            inter_center_dist = center_dists[mask].mean()
            inter_center_dist_sq = inter_center_dist
        else:
            inter_center_dist_sq = 1

        if count > 0:
            latent_dist_loss = latent_dist_loss / count / inter_center_dist_sq
        else:
            latent_dist_loss = 0.0

        loss = loss_recon + alpha * latent_dist_loss
        loss.backward()
        optimizer.step()
        if (epoch + 1) % 10 == 0 or epoch == 0:
            print(f"Epoch [{epoch+1}/{num_epochs}], Recon Loss: {loss_recon.item():.6f}, LatentDist Loss: {latent_dist_loss:.6f}, Total Loss: {loss.item():.6f}")

    encoder.eval()
    decoder.eval()
    with torch.no_grad():
        Xlatent = encoder(Xdata_train)
        Ylatent = encoder(Ydata)


    W = torch.randn(Xdata.shape[1], Ydata.shape[1]).to(device)
    W = Variable(W, requires_grad=True)

    optimizer = torch.optim.Adam([W], lr=0.01)
    lastloss = 999999

    for i in range(1000000000):
        optimizer.zero_grad()
        Wp = f.softmax(W, dim=1)
        Ypre = torch.mm(Xlatent, Wp)
        CCClose = 0
        if D != 0:
            L = torch.mm(torch.mm(Wp.t(), Idata), Wp)
            L = (L - torch.mean(L)) / torch.std(L) - 2
            L = f.sigmoid(L / 2)
            CCClose = f.cosine_similarity(Ddata, L, dim=0).mean() * 100

        Similarity1 = f.cosine_similarity(Ylatent, Ypre, dim=0).mean() * 100
        Similarity2 = f.cosine_similarity(Ylatent, Ypre, dim=1).mean() * 100
        ReguW = (
            torch.sum(Wp * torch.log(Wp), dim=1).mean() / np.log(Xdata.shape[1]) + 1
        ) * 100
        loss = -0.5 * Similarity1 - 0.5 * Similarity2 - 0.01 * ReguW - 0.004 * CCClose
        loss.backward()
        optimizer.step()
        if i % 100 == 0:
            if (lastloss - loss) / abs(lastloss) < 0.001:
                break
            else:
                lastloss = loss

                print(i)
                print(loss)
                print(Similarity1 / 100)
                print(Similarity2 / 100)
                print(CCClose / 100)
                print(ReguW / 100)

    Wp_np = Wp.cpu().detach().numpy()
    Wp_np_df = pd.DataFrame(Wp_np, index=scName_part, columns=spName_part)
    Wp_np_df_full = pd.DataFrame(0, index=scName_full, columns=spName_full).astype(float)
    Wp_np_df_full.loc[scName_part, spName_part] = Wp_np_df
    Wp_np_df_full = round(Wp_np_df_full * 1000) / 1000
    if len(genelist) > 0:
        Wp_np_df_full.to_csv(outpath + "/W.csv")
    else:
        Wp_np_df_full.to_csv(outpath + "/W_FullGene.csv")

    ### Predict Expr

    W = Wp_np_df_full
    PY = X_ori @ W

    UnionGenes = list(set(list(set(Xallgenes) | set(Yallgenes))))
    commonGenes = list(set(list(set(Xallgenes) & set(Yallgenes))))
    colsumY = np.sum(Y_ori.loc[commonGenes,], axis=0)
    colsumPY = np.sum(PY.loc[commonGenes,], axis=0)
    ratio = colsumY / (colsumPY + 0.0001)

    colsumW = np.sum(W, axis=0)
    colsumW[colsumW < 0.1] = 0
    colsumW[colsumW >= 0.1] = 1
    badcol = colsumW[colsumW == 0].index.tolist()

    Ybig = pd.DataFrame(
        0,
        index=UnionGenes,
        columns=spName_full,
    ).astype(float)
    Ytemp = Ybig.copy()
    Ytemp.loc[Yallgenes, spName_full] = Y_ori
    Ybig.loc[Xallgenes, spName_part] = PY * np.array(ratio)
    # Ybig[Ytemp > 3] = Ytemp[Ytemp > 3]
    Ybig.loc[:, badcol] = Ytemp.loc[:, badcol]
    Ybig = round(Ybig * 100) / 100
    if len(genelist) > 0:
        Ybig.to_csv(outpath + "/Ypredicted.csv")
    else:
        Ybig.to_csv(outpath + "/Ypredicted_FullGene.csv")
