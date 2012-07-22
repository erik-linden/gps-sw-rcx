function codehist = CACODEGN(PRNNo)

%  Generate the G1 and G2 sequences.
G1hist = ones(1,10);
G2hist = ones(1,10);
codehist = zeros(1023,1);
%pick appropriate taps for particular PRN
switch PRNNo
    case 1
        g2tap1 = 2; g2tap2 = 6;
    case 2
        g2tap1 = 3; g2tap2 = 7;
    case 3
        g2tap1 = 4; g2tap2 = 8;
    case 4
        g2tap1 = 5; g2tap2 = 9;
    case 5
        g2tap1 = 1; g2tap2 = 9;
    case 6
        g2tap1 = 2; g2tap2 = 10;
    case 7
        g2tap1 = 1; g2tap2 = 8;
    case 8
        g2tap1 = 2; g2tap2 = 9;
    case 9
        g2tap1 = 3; g2tap2 = 10;
    case 10
        g2tap1 = 2; g2tap2 = 3;
    case 11
        g2tap1 = 3; g2tap2 = 4;
    case 12
        g2tap1 = 5; g2tap2 = 6;
    case 13
        g2tap1 = 6; g2tap2 = 7;
    case 14
        g2tap1 = 7; g2tap2 = 8;
    case 15
        g2tap1 = 8; g2tap2 = 9;
    case 16
        g2tap1 = 9; g2tap2 = 10;
    case 17
        g2tap1 = 1; g2tap2 = 4;
    case 18
        g2tap1 = 2; g2tap2 = 5;
    case 19
        g2tap1 = 3; g2tap2 = 6;
    case 20
        g2tap1 = 4; g2tap2 = 7;
    case 21
        g2tap1 = 5; g2tap2 = 8;
    case 22
        g2tap1 = 6; g2tap2 = 9;
    case 23
        g2tap1 = 1; g2tap2 = 3;
    case 24
        g2tap1 = 4; g2tap2 = 6;
    case 25
        g2tap1 = 5; g2tap2 = 7;
    case 26
        g2tap1 = 6; g2tap2 = 8;
    case 27
        g2tap1 = 7; g2tap2 = 9;
    case 28
        g2tap1 = 8; g2tap2 = 10;
    case 29
        g2tap1 = 1; g2tap2 = 6;
    case 30
        g2tap1 = 2; g2tap2 = 7;
    case 31
        g2tap1 = 3; g2tap2 = 8;
    case 32
        g2tap1 = 4; g2tap2 = 9;
end
%  Load all rows and create C/A code
for jj = 1:1023
    %create code for step jj first
    %note that modulo 2 addition is equivalent to xor
    codehist(jj) = mod(G1hist(10)+G2hist(g2tap1)+G2hist(g2tap2),2);
    %  Create next row of G1
    G10 = xor(G1hist(3),G1hist(10));
    %  Insert next row of G1
    G1hist(:) = [G10,G1hist(1:9)];
    %  Create next row of G2
    G20 = mod(G2hist(2)+G2hist(3)+G2hist(6)+...
        G2hist(8)+G2hist(9)+G2hist(10), 2);
    %  Insert next row of G2
    G2hist(:) = [G20,G2hist(1:9)];
end


