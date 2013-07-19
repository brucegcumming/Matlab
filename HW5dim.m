%%%% Problem 3 plotting Ds rather than gammas and lambdas
nreps=5;
dims=[5];
for dim=dims;
    for rep=1:nreps
        DATAtrain = makesamples(10^3,dim);
        Xtrain=[DATAtrain.X1; DATAtrain.X2];
        Ytrain=[ones(DATAtrain.nnoise,1).*(-1) ; ones(DATAtrain.nsignal,1)]; 
        counter=0;
        DATAtest = makesamples(10^4,dim);     %generate test data
        Xtest=[DATAtest.X1;DATAtest.X2];
        Ytest=[ones(DATAtest.nnoise,1).*(-1); ones(DATAtest.nsignal,1)]; 
        gammas =[10^-4 10^-2 1 10^2];
        lambdas=[10^-7 10^-3 1 10^3];
        for i=2
            for j=4
                gamma=gammas(i);
                lambda=lambdas(j);
                S = svmtrain(Ytrain,Xtrain,['-g ' num2str(gamma) ' -c ' num2str(lambda)]); % fit the SVM
                W = S.sv_coef*S.Label(1); % these are the nonzero weights times the class labels, Mx1 vector
                SV = full(S.SVs); % these are the support vectors (training data with nonzero weights), Mxd matrix
                x=0;
                distance=0;
                for d=1:dim
                    D{d}=bsxfun(@minus,Xtest(:,d),SV(:,d)');
                    D{d}=D{d}.^2;
                    distance=distance+D{d};
                end
                K=exp(-gamma.*(distance));
                R=bsxfun(@times,K,W');
                H=sum(R,2);
                [hs x] = sort(H);
                for k = 1:length(hs)
                    upto=x(1:k);
                    ups=length(find(Ytest(upto)==1));
                    downs=length(find(Ytest(upto)==-1));
                    DR(dim,rep,k)= (DATAtest.nsignal - ups) / DATAtest.nsignal;
                    FAR(dim,rep,k)=(DATAtest.nnoise - downs) / DATAtest.nnoise;
                end
            end
        end
    end
end
figure
count=0;
for i=1:length(dims)
    dim=dims(i);
    count=count+1;
    subplot(length(dims),1,count);
    plot(squeeze(FAR(i,:,:))',squeeze(DR(i,:,:))','HandleVisibility','off');
    title(['# dimension = ',num2str(dim)]);     
    set(get(gca,'XLabel'),'string','False Alarm Rate');
    set(get(gca,'YLabel'),'string','Detection Rate');
end 