% Rajat Kumar
% 2021185
t_initial = 0;

% no. of paths
K = 4;

%number of sinusoids used to model a single path.
Nk = 48;

% maximum doppler shift
fd = 0;

% Limits
N1=0;
N2=15;

% Sampling Frequencies
sampling_freq_Hz=[1e6,10e3,1e3];

%Sampling Time
Ts=zeros(1,length(sampling_freq_Hz));
for i=1:length(sampling_freq_Hz)
    Ts(i)=1/sampling_freq_Hz(i);
end

% Path Delays
Tau_s=[0,5e-6,10e-6,15e-6];

% 2D frequencies array to store in phase component and quadrature
% components
frequencies = zeros(2,K);

% We Know that each element of frequencies array will be computed by the
% formula f(K,n) = fd*cos((pi/2Nk)*(n-1/2)+aK) but fd=0 in our
% case. Therefore frequencies array will contain zeroes.

% Constellation points given in assignment
constellation = [0.7+0.7i; -0.7+0.7i; -0.7-0.7i; 0.7-0.7i];

% Plotting the required parts iterating over each sampling frequency
for j = 1:length(Ts)
tau_s=Tau_s/Ts(j);

% randomly generating theta
theta = zeros(2, K, Nk);
val=0;
for i=1:K
    for l=1:Nk
        val=2*pi*rand(1);
        theta(1,i,l)=val;
        theta(2,i,l)=val+pi/2;
    end
end

% mu(K) = sqrt(2/Nk)*summation(cos(2*pi*fK,n*t + theta K,n)). But fK,n
% is zero in our case. therefore,  mu(K) = sqrt(2/Nk)*summation(cos(theta K,n))
mu = zeros(2, K);
for i = 1:2
    for l = 1:K
        sum=0;
        for x=1:Nk
            sum=sum+ cos(theta(i,l,x));
        end
        sum= sqrt(2/Nk)*sum;
        mu(i,l)=sum;
    end
end

% Computing fading waveform for each path (zk)
zk = mu(1,:) + 1j*mu(2,:);

% we know that ak=sqrt(Omega_k)*zk
% we can assume Omega_k=1
% therefore ak=zk
% For Omega_k=1 E[|ak|^2]=1 and we know that E[constant]=constant
% Which implies |ak|=|zk|=1

% ak=zk according to the correction made in assignment 3
ak=zeros(1,K);
for i=1:K
%     ak(i)=zk(i)/abs(zk(i));
    ak(i)=zk(i);
end

% Initialize gn array
gn = zeros(1, N2+1); 

for n = 1:N2+1

    % Initialize sum for the current time index n
    sum = 0;
    
    % Loop over each path k to compute the contribution to gn(n)
    for k = 1:K
        % Calculate the argument of the sinc function
        arg_sinc = tau_s(k) - (n - 1);
        
        % Evaluate the sinc function for the current path k
        sinc_value = my_sinc(arg_sinc);
        
        % Compute the contribution of the current path to gn(n)
        contribution_k = ak(k) * sinc_value;
        
        % Accumulate the contribution to the sum
        sum = sum + contribution_k;
    end
    
    % Assign the computed sum to the n-th element of gn
    gn(n) = sum;
end

    % Generate 1024 QPSK symbols
    num_symbols = 1024;
    qpsk_symbols = constellation(randi([1, 4], num_symbols, 1));
    
    % Modify the figure size
    figure('Position', [100, 100, 1400, 1000]); 
    
    % Plotting Constellation points and generated QPSK symbols (Part 2)
    subplot(3,2,1);
    scatter(real(qpsk_symbols), imag(qpsk_symbols), 'b', 'DisplayName', 'Generated Symbols');
    hold on;
    scatter(real(constellation), imag(constellation), 'g', 'filled', 'DisplayName', 'Constellation Points');
    hold off;
    title(['Scatter Plot of QPSK Symbols (Fs = ', num2str(sampling_freq_Hz(j)), ' Hz)']);
    xlabel('In-phase');
    ylabel('Quadrature');
    axis square;
    grid on;
    legend;

    % Pass QPSK symbols through the channel
    received_symbols = conv(qpsk_symbols, gn);
    received_symbols=received_symbols(1:1024);

    % Plot scatter plot of received symbols overlaid with transmitted
    % symbols (Part 3)
    subplot(3,2,2);
    scatter(real(received_symbols), imag(received_symbols), 'r',  'DisplayName', 'Received Symbols');
    hold on;
    scatter(real(constellation), imag(constellation), 'g', 'filled', 'DisplayName', 'Constellation Points');
    scatter(real(qpsk_symbols), imag(qpsk_symbols), 'b', 'DisplayName', 'Transmitted Symbols');
    hold off;
    title(['Received and Transmitted QPSK Symbols (Fs = ', num2str(sampling_freq_Hz(j)), ' Hz)']);
    xlabel('In-phase');
    ylabel('Quadrature');
    axis square;
    grid on;
    legend;

    % Divide the received symbols by the first channel coefficient
    received_symbols = received_symbols / gn(1);
    
    % Plotting received signals after dividing by gn(1) (Part 4)
    subplot(3,2,3);
    scatter(real(received_symbols), imag(received_symbols), 'r',  'DisplayName', 'Received Symbols');
    hold on;
    scatter(real(constellation), imag(constellation), 'g', 'filled', 'DisplayName', 'Constellation Points');
    scatter(real(qpsk_symbols), imag(qpsk_symbols), 'b',  'DisplayName', 'Transmitted Symbols');
    hold off;
    title(['Received and Transmitted QPSK Symbols after dividing by first channel coeff (Fs = ', num2str(sampling_freq_Hz(j)), ' Hz)']);
    xlabel('In-phase');
    ylabel('Quadrature');
    axis square;
    grid on;
    legend;
    
    % OFDM
    
    % Dividing into parallel streams
    Parallel_qpsk = cell(1, 16);
    parall=[];
    y=0;
    for i=1:16
        parall=zeros(1,64);
        for x=1:64
            y=y+1;
            parall(x)=qpsk_symbols(y);
        end
        Parallel_qpsk{i}=parall;
    end
    
    % performing IFFT
    for i=1:16
        temp=Parallel_qpsk{i};
        Parallel_qpsk{i}=ifft(temp);
    end
    
    % Adding Cyclic Coefficients
    for i=1:16
        parall=Parallel_qpsk{i};
        cyclic_coeff=parall(50:64);
        new_parall=[cyclic_coeff,parall];
        Parallel_qpsk{i}=new_parall;
    end
    
    % Parallel to serial conversion
    Serial=[];
    for i=1:16
        Serial=[Serial,Parallel_qpsk{i}];
    end
    
    % Passing through channel
    Serial_output=conv(Serial,gn);
    Serial=Serial_output(1:1264);
    
    % Serial to Parallel
    y=0;
    for i=1:16
        parall=zeros(1,79);
        for x=1:79
            y=y+1;
            parall(x)=Serial(y);
        end
        Parallel_qpsk{i}=parall;
    end
    
    % Removing Cyclic coefficients from each of these parallel streams
    for i =1:16
        temp=Parallel_qpsk{i}(16:79);
        Parallel_qpsk{i}=temp;
    end
    
    % Performing FFT
    for i=1:16
        temp=Parallel_qpsk{i};
        Parallel_qpsk{i}=fft(temp);
    end
    
    % Parallel to serial conversion
    Ofdm_output=[];
    for i=1:16
        Ofdm_output=[Ofdm_output,Parallel_qpsk{i}];
    end
    
    % Plotting received OFDM symbols (Part 6)
    subplot(3,2,4);
    scatter(real(Ofdm_output), imag(Ofdm_output), 'r',  'DisplayName', 'Received Symbols');
    hold on;
    scatter(real(constellation), imag(constellation), 'g', 'filled', 'DisplayName', 'Constellation Points');
    scatter(real(qpsk_symbols), imag(qpsk_symbols), 'b',  'DisplayName', 'Transmitted Symbols');
    hold off;
    title('Transmitted QPSK Symbols and demodulated OFDM Symbols');
    xlabel('In-phase');
    ylabel('Quadrature');
    axis square;
    grid on;
    legend;
    
    % Perform N-point FFT of the channel
    Hf = fft(gn, 64);
        
    % Dividing ofdm_output by Hf (element-wise)
    x = 1;

    % Final OFDM output 
    ofdm_serial = [];
    for i = 1:16
        z = x + 63;
        par = Ofdm_output(x:z);
        % Add a check for zero values in Hf
        if any(Hf == 0)
            % Handle the case where Hf contains zero values
            Hf_nonzero = Hf + (Hf == 0) * eps;
            par = par ./ Hf_nonzero;
        else
            par = par ./ Hf;
        end
        x = z + 1;
        ofdm_serial = [ofdm_serial, par];
    end
    
    
     % Scatter plot of demodulated symbols divided by Hf (Part 7)
        subplot(3,2,5);
        scatter(real(ofdm_serial), imag(ofdm_serial), 'r',  'DisplayName', 'Received Symbols (After Division)');
        hold on;
        scatter(real(constellation), imag(constellation), 'g', 'filled', 'DisplayName', 'Constellation Points');
        scatter(real(qpsk_symbols), imag(qpsk_symbols), 'b', 'DisplayName', 'Transmitted Symbols');
        hold off;
        title(['Received Symbols (After Dividing OFDM by Hf) vs Transmitted Symbols (Fs = ', num2str(sampling_freq_Hz(j)), ' Hz)']);
        xlabel('In-phase');
        ylabel('Quadrature');
        axis square;
        grid on;
        legend;
end
% Define sinc function
function y = my_sinc(x)
    y = sin(pi * x) ./ (pi * x);
    y(x == 0) = 1;  % Set the value at x=0 to 1 manually to avoid division by zero
end