
## Log into Event Engine

1. Go to https://dashboard.eventengine.run.

2. Enter the event hash (provided at the beginning of this workshop). Then click **Accept Terms and Login**
![](images/eventengine1.png)

3. Select **Sign In with your social account**
![](images/eventengine2.png)

4. Enter your amazon.com credentials. These credentials are only used to identify you for this workshop.
![](images/eventengine3.png)

5. Once the event has started, the **Team Dashboard** will become available.
![](images/eventengine4.png)

6. Click on **SSH Key**, then **Download Key**. Save this key somewhere safe.
![](images/eventengine5.png)

7. Click **OK** to close the window and then select **AWS Console** from the main Dashboard.
![](images/eventengine6.png)

8. Click **Open AWS Console** to open the AWS Console
![](images/eventengine7.png)



## Launch the CloudFormation Stack

[![cloudformation-launch-stack](images/launchstack.png)](https://console.aws.amazon.com/cloudformation/home?region=us-west-2#/stacks/new?stackName=RStudio&templateURL=https://uc-aws-workshops.s3-us-west-2.amazonaws.com/rstudio-sslv5.yaml)

### Instructions
1. Launch the AWS CloudFormation in your AWS account using the **Launch Stack** button above.
![](images/cf_create_stack.png)
2. Click "Next"
3. Enter Parameters for the stack.
    - For the **KeyPair**, select "ee-default-keypair" from the dropdown
    - For the **UserList**, choose a username and password, and enter them in the format "username,password"
    - For the **VPCId** parameter, use your Default VPC (172.31.0.0/16)
    - For the **VPCSubnet** parameter, choose a subnet within the Default VPC (172.31.0.0/20)
    - Leave the other defaults alone.
![](images/cf_stack_parameters.png)

4. Click **Next**
5. Click **Next**
6. Acknowledge that AWS CloudFormation might create IAM resources, and click **Create Stack**
![](images/cf_launch_stack.png)

7. After 5 minutes, follow the link in the **Outputs** tab of your AWS CloudForamtion Stack to access RStudio.
![](images/cf_resources.png)

8. Accept the warning from your browser about the certificate being self-signed.  This gives us encrypted, HTTPS access to RStudio without purchasing a domain name or SSL certificate.

9. Login to RStudio using the credentials you provided to the AWS CloudFormation template.
!()[images/rstudio1.png]
